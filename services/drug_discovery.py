import hashlib
import random
from collections import Counter
from pathlib import Path

from Bio.Seq import Seq


def _sanitize_dna(sequence):
    cleaned = "".join(ch for ch in str(sequence or "").upper() if ch in {"A", "C", "G", "T", "N"})
    return cleaned.replace("N", "")


def _gc_content(dna_sequence):
    if not dna_sequence:
        return 0.0
    gc = dna_sequence.count("G") + dna_sequence.count("C")
    return round((gc / len(dna_sequence)) * 100, 2)


def _reverse_complement(dna_sequence):
    return str(Seq(dna_sequence).reverse_complement())


def _codon_usage_summary(dna_sequence, top_n=20):
    codons = [dna_sequence[i:i + 3] for i in range(0, len(dna_sequence) - 2, 3)]
    codons = [c for c in codons if len(c) == 3]
    counts = Counter(codons)
    total = sum(counts.values()) or 1
    sorted_pairs = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)[:top_n]
    summary = {
        codon: {
            "count": count,
            "frequency": round((count / total) * 100, 2),
        }
        for codon, count in sorted_pairs
    }
    return summary


def get_best_protein(dna_sequence, frame=None):
    """Find the longest translated ORF using translate(to_stop=True)."""
    sequence = _sanitize_dna(dna_sequence)
    if not sequence:
        return {
            "frame": 0,
            "protein": "",
            "protein_length": 0,
            "protein_preview": "",
            "stop_position": None,
        }

    frames = [frame] if frame in {0, 1, 2} else [0, 1, 2]
    best = {
        "frame": 0,
        "protein": "",
        "protein_length": 0,
        "protein_preview": "",
        "stop_position": None,
    }

    for frame_idx in frames:
        frame_best = {
            "frame": frame_idx,
            "protein": "",
            "protein_length": 0,
            "protein_preview": "",
            "stop_position": None,
        }

        for nt_pos in range(frame_idx, len(sequence) - 2, 3):
            if sequence[nt_pos:nt_pos + 3] != "ATG":
                continue
            sub_seq = sequence[nt_pos:]
            trimmed = sub_seq[: len(sub_seq) - (len(sub_seq) % 3)]
            if not trimmed:
                continue

            try:
                translated = str(Seq(trimmed).translate(to_stop=True))
            except Exception:
                continue

            if len(translated) > frame_best["protein_length"]:
                stop_nt = nt_pos + (len(translated) * 3) + 3
                frame_best = {
                    "frame": frame_idx,
                    "protein": translated,
                    "protein_length": len(translated),
                    "protein_preview": translated[:50],
                    "stop_position": stop_nt if stop_nt <= len(sequence) else None,
                }

        if frame_best["protein_length"] > best["protein_length"]:
            best = frame_best

    return best


def analyze_dna(sequence):
    """Run DNA analysis with ORF scan across all 3 reading frames."""
    dna_sequence = _sanitize_dna(sequence)
    if not dna_sequence:
        raise ValueError("Empty sequence after cleanup. Upload or provide valid DNA bases (A/C/G/T).")
    if len(dna_sequence) < 60:
        raise ValueError("Sequence is too short for meaningful ORF and docking simulation (minimum 60 nt).")

    gc = _gc_content(dna_sequence)
    rev = _reverse_complement(dna_sequence)

    frame_orfs = []
    for frame_idx in [0, 1, 2]:
        frame_orfs.append(get_best_protein(dna_sequence, frame=frame_idx))

    best_orf = max(frame_orfs, key=lambda item: item["protein_length"])
    aa_freq = Counter(best_orf["protein"])

    return {
        "gc_content": gc,
        "length": len(dna_sequence),
        "sequence": dna_sequence,
        "reverse_complement": rev,
        "orf_frames": frame_orfs,
        "protein_length": best_orf["protein_length"],
        "protein_preview": best_orf["protein_preview"],
        "protein_full": best_orf["protein"],
        "stop_position": best_orf["stop_position"],
        "best_frame": best_orf["frame"],
        "amino_acid_frequency": dict(sorted(aa_freq.items(), key=lambda kv: kv[1], reverse=True)),
        "codon_usage_summary": _codon_usage_summary(dna_sequence, top_n=20),
    }


def simulate_binding_energy(dna_analysis, ligand_name=None):
    """Generate multiple docking-like poses with controlled, sequence-driven variation."""
    gc_content = float(dna_analysis.get("gc_content", 0.0))
    seq_length = int(dna_analysis.get("length", 0))
    protein_length = int(dna_analysis.get("protein_length", 0))

    seed_text = (
        f"{gc_content:.2f}|{seq_length}|{protein_length}|"
        f"{ligand_name or ''}|{dna_analysis.get('protein_preview', '')}"
    )
    seed = int(hashlib.sha256(seed_text.encode("utf-8", errors="ignore")).hexdigest()[:12], 16)
    rng = random.Random(seed)

    # In real docking, different poses come from different ligand orientations/conformations.
    pose_count = rng.randint(10, 15)

    # This base value is a heuristic approximation, not a physics-based docking score.
    base_energy = -6.2 - ((gc_content - 50.0) / 40.0)
    base_energy -= min(seq_length / 15000.0, 1.2)
    base_energy -= min(protein_length / 1200.0, 1.0)
    base_energy = max(min(base_energy, -4.5), -9.2)

    energies = []
    for pose_index in range(pose_count):
        pose_shift = (pose_index - (pose_count / 2.0)) * 0.03
        local_noise = rng.uniform(-0.42, 0.42)
        energy = base_energy + pose_shift + local_noise
        energies.append(round(max(min(energy, -4.0), -10.0), 2))

    best_energy = min(energies)
    best_index = energies.index(best_energy)
    average_energy = round(sum(energies) / len(energies), 2)

    return {
        "binding_energies": energies,
        "best_energy": best_energy,
        "average_energy": average_energy,
        "best_pose": best_index + 1,
    }


def generate_graph(binding_energies):
    """Return graph-ready line-plot payload for pose vs binding energy."""
    if not binding_energies:
        return {
            "labels": [],
            "values": [],
            "best_pose": None,
            "best_energy": None,
            "x_label": "Pose Number",
            "y_label": "Binding Energy (kcal/mol)",
        }

    best_energy = min(binding_energies)
    best_pose = binding_energies.index(best_energy) + 1
    return {
        "labels": list(range(1, len(binding_energies) + 1)),
        "values": binding_energies,
        "best_pose": best_pose,
        "best_energy": best_energy,
        "x_label": "Pose Number",
        "y_label": "Binding Energy (kcal/mol)",
    }


def _build_interaction_summary(best_energy):
    # Lower binding energy typically implies stronger predicted affinity in docking conventions.
    if best_energy <= -8.5:
        return "High-affinity pose pattern with strong predicted target engagement."
    if best_energy <= -7.0:
        return "Moderate-to-strong predicted affinity across sampled poses."
    return "Weak-to-moderate predicted affinity; optimization is likely needed."


def _build_target_name(sequence_name, ligand_name):
    base = sequence_name or ligand_name or "UploadedSequence"
    stem = Path(base).stem.strip().replace("_", " ")
    return stem if stem else "UploadedSequence"


def _infer_file_type(file_name):
    if not file_name or "." not in file_name:
        return "unknown"
    return file_name.rsplit(".", 1)[-1].lower()


def _sequence_fingerprint(sequence):
    digest = hashlib.sha256(sequence.encode("utf-8", errors="ignore")).hexdigest()
    return digest[:16]


def _estimate_bond_profile(result):
    """Estimate docking interaction details from simulated affinity and sequence context."""
    best_energy = float(result.get("best_energy", -6.0))
    protein_length = int(result.get("protein_length", 0))
    seed_text = (
        f"{result.get('sequence_fingerprint', '')}|{result.get('ligand_name', '')}|"
        f"{best_energy}|{protein_length}"
    )
    seed = int(hashlib.sha256(seed_text.encode("utf-8", errors="ignore")).hexdigest()[:12], 16)
    rng = random.Random(seed)

    affinity_scale = max(0.0, min(((-best_energy) - 4.0) / 6.0, 1.0))

    h_bonds = int(round(2 + (affinity_scale * 6) + rng.uniform(-1, 1)))
    hydrophobic_contacts = int(round(6 + (affinity_scale * 14) + rng.uniform(-2, 2)))
    ionic_interactions = int(round(1 + (affinity_scale * 4) + rng.uniform(-1, 1)))
    pi_stacking = int(round(0 + (affinity_scale * 3) + rng.uniform(-1, 1)))

    h_bonds = max(0, h_bonds)
    hydrophobic_contacts = max(1, hydrophobic_contacts)
    ionic_interactions = max(0, ionic_interactions)
    pi_stacking = max(0, pi_stacking)

    pocket_start = max(8, int((protein_length * 0.18) + rng.randint(-8, 8)))
    pocket_end = max(pocket_start + 8, int((protein_length * 0.36) + rng.randint(-8, 8)))

    return {
        "hydrogen_bonds": h_bonds,
        "hydrophobic_contacts": hydrophobic_contacts,
        "ionic_interactions": ionic_interactions,
        "pi_stacking": pi_stacking,
        "predicted_binding_pocket": f"AA {pocket_start}-{pocket_end}",
    }


def _build_report(result):
    bond = result.get("bond_profile", {})
    return "\n".join([
        "========================================",
        "   DNA-GUIDED DOCKING SIMULATION REPORT ",
        "========================================",
        "",
        f"Target: {result['target_name']}",
        f"DNA File: {result['protein_name']}",
        f"DNA Fingerprint: {result['sequence_fingerprint']}",
        f"Ligand File: {result['ligand_name']} ({result['ligand_type']})",
        f"Sequence Length: {result['length']} nt",
        f"GC Content: {result['gc_content']}%",
        f"Best ORF Frame: {result['best_frame']}",
        f"Best ORF Protein Length: {result['protein_length']} aa",
        f"Best ORF Protein Preview: {result['protein_preview'] or 'N/A'}",
        f"Best ORF Stop Position: {result['stop_position']}",
        f"DNA Sequence Preview: {result['dna_sequence_preview']}",
        "",
        f"Best Binding Energy: {result['best_energy']} kcal/mol",
        f"Average Binding Energy: {result['average_energy']} kcal/mol",
        f"Number of Poses: {len(result['binding_energies'])}",
        f"Hydrogen Bonds (estimated): {bond.get('hydrogen_bonds', 'N/A')}",
        f"Hydrophobic Contacts (estimated): {bond.get('hydrophobic_contacts', 'N/A')}",
        f"Ionic Interactions (estimated): {bond.get('ionic_interactions', 'N/A')}",
        f"Pi-stacking Events (estimated): {bond.get('pi_stacking', 'N/A')}",
        f"Predicted Binding Pocket: {bond.get('predicted_binding_pocket', 'N/A')}",
        "",
        "Scientific Note:",
        "Multiple poses are sampled because ligands can dock in different orientations.",
        "Lower (more negative) binding energy is interpreted as stronger predicted affinity.",
        "This is a simulation-based approximation of docking behavior, not a replacement for AutoDock/Vina or wet-lab assays.",
    ])


def simulate_docking(
    disease=None,
    has_protein=False,
    has_ligand=False,
    protein_name=None,
    ligand_name=None,
    protein_sequence=None,
    force_uploaded_mode=False,
):
    """Main orchestrator that combines DNA analysis and multi-pose affinity simulation."""
    _ = disease, has_protein, has_ligand, force_uploaded_mode  # kept for backward-compatible calls
    dna_analysis = analyze_dna(protein_sequence)
    binding = simulate_binding_energy(dna_analysis, ligand_name=ligand_name)
    graph = generate_graph(binding["binding_energies"])

    result = {
        "target_name": _build_target_name(protein_name, ligand_name),
        "gc_content": dna_analysis["gc_content"],
        "length": dna_analysis["length"],
        "reverse_complement": dna_analysis["reverse_complement"],
        "orf_frames": dna_analysis["orf_frames"],
        "protein_length": dna_analysis["protein_length"],
        "best_frame": dna_analysis["best_frame"],
        "protein_preview": dna_analysis["protein_preview"],
        "protein_full": dna_analysis["protein_full"],
        "stop_position": dna_analysis["stop_position"],
        "dna_sequence_preview": dna_analysis["sequence"][:120],
        "dna_sequence_full": dna_analysis["sequence"],
        "reverse_complement_preview": dna_analysis["reverse_complement"][:120],
        "amino_acid_frequency": dna_analysis["amino_acid_frequency"],
        "codon_usage_summary": dna_analysis["codon_usage_summary"],
        "binding_energies": binding["binding_energies"],
        "pose_count": len(binding["binding_energies"]),
        "best_energy": binding["best_energy"],
        "average_energy": binding["average_energy"],
        "best_pose": binding["best_pose"],
        "graph": graph,
        "interaction": _build_interaction_summary(binding["best_energy"]),
        "lower_is_better": True,
        "protein_name": protein_name or "Uploaded DNA",
        "ligand_name": ligand_name or "Uploaded Ligand",
        "ligand_type": _infer_file_type(ligand_name),
        "sequence_fingerprint": _sequence_fingerprint(dna_analysis["sequence"]),
        "sequence_length": dna_analysis["length"],
        # Backward compatibility keys used by existing UI/paper builder.
        "binding_energy": binding["best_energy"],
        "disease": _build_target_name(protein_name, ligand_name),
        "demo_mode": False,
    }

    result["bond_profile"] = _estimate_bond_profile(result)

    result["report"] = _build_report(result)
    return result
