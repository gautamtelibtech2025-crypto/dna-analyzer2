"""
app.py - DNA & Protein Bioinformatics Analyzer
Flask backend that fetches DNA/Protein sequences from NCBI and performs analysis using BioPython.

Routes:
  GET  /                      → Serve the main page
  POST /api/search             → Search NCBI nucleotide DB for a gene name
  POST /api/analyze            → Fetch DNA sequence and return full analysis
  POST /api/search-protein     → Search NCBI protein DB for a protein name
  POST /api/analyze-protein    → Fetch protein sequence and return full analysis
  POST /api/download-pdf       → Generate and download a full PDF report
"""

from flask import Flask, render_template, request, jsonify, send_file
from Bio import Entrez, SeqIO, Medline
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import Counter
import re
import io
import json
import os
import socket
import requests
from dotenv import load_dotenv
import time
from services.drug_discovery import analyze_dna as analyze_dna_core, simulate_docking
app = Flask(__name__)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
load_dotenv(os.path.join(BASE_DIR, ".env"))

# NCBI requires a REAL email for Entrez API usage
Entrez.email = os.getenv("ENTREZ_EMAIL", "gautamrathore177@gmail.com")
GEMINI_API_KEY = os.getenv("GEMINI_API_KEY")
GEMINI_MODEL = os.getenv("GEMINI_MODEL", "gemini-2.0-flash")


def _gemini_model_candidates():
    models = [
        GEMINI_MODEL,
        "gemini-2.5-flash",
        "gemini-2.5-flash-lite",
        "gemini-2.5-pro",
        "gemini-2.0-flash",
        "gemini-2.0-flash-001",
        "gemini-2.0-flash-lite-001",
        "gemini-1.5-flash-latest",
        "gemini-1.5-flash",
    ]
    seen = set()
    ordered = []
    for m in models:
        if not m or m in seen:
            continue
        seen.add(m)
        ordered.append(m)
    return ordered

# Set a global socket timeout so NCBI Entrez calls never hang forever (30s)
socket.setdefaulttimeout(30)

# Simple in-memory cache to avoid re-generating the same topic repeatedly.
RESEARCH_CACHE_TTL_SECONDS = 60 * 60
research_cache = {}
SIM_ENRICH_CACHE_TTL_SECONDS = 60 * 60
simulation_enrichment_cache = {}

DEMO_PROJECTS = {
    "brca1_demo": {
        "dna_name": "BRCA1_demo.fasta",
        "ligand_name": "olaparib_demo.mol",
        "fasta_header": ">BRCA1 Homo sapiens demo sequence",
        "dna_sequence": "ATG" + ("GCC" * 120) + "TAA" + ("GGT" * 50),
    },
    "tp53_demo": {
        "dna_name": "TP53_demo.fasta",
        "ligand_name": "nutlin3_demo.mol",
        "fasta_header": ">TP53 Homo sapiens demo sequence",
        "dna_sequence": "ATG" + ("GAA" * 110) + "TGA" + ("CCT" * 45),
    },
    "egfr_demo": {
        "dna_name": "EGFR_demo.fasta",
        "ligand_name": "gefitinib_demo.mol",
        "fasta_header": ">EGFR Homo sapiens demo sequence",
        "dna_sequence": "ATG" + ("ACC" * 130) + "TAG" + ("GAT" * 40),
    },
}


# ─────────────────────────────────────────────
# Bioinformatics Analysis Functions
# ─────────────────────────────────────────────

def search_gene(gene_name, max_results=3):
    """
    Search the NCBI Nucleotide database for a gene name.
    Returns a list of dicts with id, title, length for the top results.
    """
    try:
        # Search specifically in Homo sapiens (human) for relevant results
        handle = Entrez.esearch(
            db="nucleotide",
            term=f"{gene_name}[Gene Name] Homo sapiens[Organism]",
            retmax=max_results
        )
        record = Entrez.read(handle)
        handle.close()

        ids = record.get("IdList", [])
        if not ids:
            return []

        # Fetch summaries (title, length, etc.) for each hit
        summary_handle = Entrez.esummary(db="nucleotide", id=",".join(ids))
        summaries = Entrez.read(summary_handle)
        summary_handle.close()

        results = []
        for s in summaries:
            results.append({
                "id": s["Id"],
                "title": s["Title"],
                "length": s.get("Length", "?"),
            })
        return results

    except Exception as e:
        raise RuntimeError(f"NCBI search failed: {str(e)}")


def fetch_sequence(ncbi_id):
    """
    Fetch the full GenBank record from NCBI by its accession/ID.
    Returns a BioPython SeqRecord object.
    """
    try:
        handle = Entrez.efetch(
            db="nucleotide", id=ncbi_id,
            rettype="gb", retmode="text"
        )
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return record
    except Exception as e:
        raise RuntimeError(f"Failed to fetch sequence: {str(e)}")


def compute_gc_content(seq):
    """Returns GC content as a percentage, e.g. 52.34."""
    return round(gc_fraction(seq) * 100, 2)


def get_reverse_complement(seq):
    """Returns the reverse complement of the DNA sequence as a string."""
    return str(seq.reverse_complement())


def translate_sequence(seq):
    """
    Translates the DNA sequence into a protein (amino acid) sequence.
    Trims to a multiple of 3 to avoid partial-codon errors.
    Stops at the first stop codon.
    """
    trimmed = seq[: len(seq) - (len(seq) % 3)]
    protein = trimmed.translate(to_stop=True)
    return str(protein)


def codon_frequency(seq):
    """
    Count the frequency of each codon (3-letter DNA triplet).
    Returns a dict of codon → count, sorted descending, top 20.
    """
    codons = [
        str(seq[i : i + 3])
        for i in range(0, len(seq) - 2, 3)
        if len(seq[i : i + 3]) == 3
    ]
    freq = Counter(codons)
    # Sort by count descending and keep the top 20 codons
    return dict(sorted(freq.items(), key=lambda x: x[1], reverse=True)[:20])


def nucleotide_counts(seq):
    """Return individual base counts for A, T, G, C."""
    seq_str = str(seq).upper()
    return {
        "A": seq_str.count("A"),
        "T": seq_str.count("T"),
        "G": seq_str.count("G"),
        "C": seq_str.count("C"),
    }


def analyze_sequence(record):
    """
    Run all analyses on a SeqRecord and return a results dictionary.
    Includes both truncated previews (for fast UI) and full sequences.
    """
    seq = record.seq
    dna_str = str(seq)
    core = analyze_dna_core(dna_str)

    codon_frequency_counts = {
        codon: values.get("count", 0)
        for codon, values in core.get("codon_usage_summary", {}).items()
    }

    return {
        "accession": record.id,
        "description": record.description,
        "sequence_length": core["length"],
        "sequence_preview": core["sequence"][:200],
        "sequence_full": core["sequence"],
        "gc_content": core["gc_content"],
        "reverse_complement": core["reverse_complement"][:200],
        "reverse_complement_full": core["reverse_complement"],
        "protein": core["protein_preview"],
        "protein_full": core["protein_full"],
        "protein_length": core["protein_length"],
        "protein_preview": core["protein_preview"],
        "stop_position": core["stop_position"],
        "orf_frames": core["orf_frames"],
        "codon_frequency": codon_frequency_counts,
        "codon_usage_summary": core["codon_usage_summary"],
        "amino_acid_frequency": core["amino_acid_frequency"],
        "nucleotide_counts": nucleotide_counts(core["sequence"]),
        # Required flat keys for simple consumers.
        "length": core["length"],
    }


# ─────────────────────────────────────────────
# Protein Analysis Functions
# ─────────────────────────────────────────────

def search_protein(protein_name, max_results=3):
    """
    Search the NCBI Protein database for a protein name.
    Returns a list of dicts with id, title, length for the top results.
    """
    try:
        # First try strict Protein Name search
        handle = Entrez.esearch(
            db="protein",
            term=f"{protein_name}[Protein Name] AND Homo sapiens[Organism]",
            retmax=max_results
        )
        record = Entrez.read(handle)
        handle.close()

        ids = record.get("IdList", [])
        
        # Fallback to general search if strict search missed
        if not ids:
            handle = Entrez.esearch(
                db="protein",
                term=f"{protein_name} AND Homo sapiens[Organism]",
                retmax=max_results
            )
            record = Entrez.read(handle)
            handle.close()
            ids = record.get("IdList", [])
            
        if not ids:
            return []

        summary_handle = Entrez.esummary(db="protein", id=",".join(ids))
        summaries = Entrez.read(summary_handle)
        summary_handle.close()

        results = []
        for s in summaries:
            results.append({
                "id": s["Id"],
                "title": s["Title"],
                "length": s.get("Length", "?"),
            })
        return results

    except Exception as e:
        raise RuntimeError(f"NCBI protein search failed: {str(e)}")


def fetch_protein_sequence(ncbi_id):
    """
    Fetch the full protein record from NCBI by its accession/ID.
    Returns a BioPython SeqRecord object.
    """
    try:
        handle = Entrez.efetch(
            db="protein", id=ncbi_id,
            rettype="gp", retmode="text"
        )
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return record
    except Exception as e:
        raise RuntimeError(f"Failed to fetch protein sequence: {str(e)}")


def analyze_protein_sequence(record):
    """
    Run comprehensive protein analysis on a SeqRecord.
    Uses ProtParam for physicochemical properties.
    """
    seq_str = str(record.seq)

    # Remove any stop codon '*' or unknown 'X' for ProtParam
    clean_seq = seq_str.replace('*', '').replace('X', '').replace('x', '')

    try:
        analysis = ProteinAnalysis(clean_seq)
    except Exception as e:
        raise RuntimeError(f"Protein analysis failed: {str(e)}")

    # Amino acid counts
    aa_count = analysis.count_amino_acids()

    # Amino acid percentage
    aa_percent = {}
    try:
        aa_percent = {k: round(v * 100, 2) for k, v in analysis.get_amino_acids_percent().items()}
    except Exception:
        pass

    # Molecular weight
    try:
        mol_weight = round(analysis.molecular_weight(), 2)
    except Exception:
        mol_weight = "N/A"

    # Isoelectric point
    try:
        iso_point = round(analysis.isoelectric_point(), 2)
    except Exception:
        iso_point = "N/A"

    # Aromaticity
    try:
        aromaticity = round(analysis.aromaticity(), 4)
    except Exception:
        aromaticity = "N/A"

    # Instability index
    try:
        instability = round(analysis.instability_index(), 2)
        is_stable = instability < 40
    except Exception:
        instability = "N/A"
        is_stable = None

    # GRAVY (Grand Average of Hydropathicity)
    try:
        gravy = round(analysis.gravy(), 4)
    except Exception:
        gravy = "N/A"

    # Secondary structure fraction
    try:
        helix, turn, sheet = analysis.secondary_structure_fraction()
        sec_structure = {
            "helix": round(helix * 100, 2),
            "turn": round(turn * 100, 2),
            "sheet": round(sheet * 100, 2),
        }
    except Exception:
        sec_structure = {"helix": 0, "turn": 0, "sheet": 0}

    # Top 20 amino acids by count
    aa_frequency = dict(
        sorted(aa_count.items(), key=lambda x: x[1], reverse=True)[:20]
    )

    # Charge at pH 7
    try:
        charge_at_7 = round(analysis.charge_at_pH(7.0), 2)
    except Exception:
        charge_at_7 = "N/A"

    return {
        "accession": record.id,
        "description": record.description,
        "sequence_length": len(seq_str),
        "sequence_preview": seq_str[:200],
        "sequence_full": seq_str,
        "molecular_weight": mol_weight,
        "isoelectric_point": iso_point,
        "aromaticity": aromaticity,
        "instability_index": instability,
        "is_stable": is_stable,
        "gravy": gravy,
        "charge_at_ph7": charge_at_7,
        "secondary_structure": sec_structure,
        "amino_acid_count": aa_count,
        "amino_acid_percent": aa_percent,
        "amino_acid_frequency": aa_frequency,
    }


def search_pubmed_papers(topic, max_results=3):
    """
    Fetch latest PubMed papers for a topic and return core bibliographic metadata.
    """
    try:
        handle = Entrez.esearch(
            db="pubmed",
            term=f"{topic}[Title/Abstract]",
            sort="pub date",
            retmax=max_results,
        )
        search_record = Entrez.read(handle)
        handle.close()

        ids = search_record.get("IdList", [])
        if not ids:
            return []

        fetch_handle = Entrez.efetch(
            db="pubmed",
            id=",".join(ids),
            rettype="medline",
            retmode="text",
        )
        records = list(Medline.parse(fetch_handle))
        fetch_handle.close()

        papers = []
        for rec in records:
            title = str(rec.get("TI", "")).strip()
            abstract = " ".join(str(rec.get("AB", "")).split())
            dp = str(rec.get("DP", ""))
            year_match = re.search(r"\b(19|20)\d{2}\b", dp)
            year = year_match.group(0) if year_match else "N/A"
            pmid = str(rec.get("PMID", "")).strip()
            journal = str(rec.get("JT", "")).strip() or str(rec.get("TA", "")).strip()
            authors = [str(a).strip() for a in rec.get("AU", []) if str(a).strip()]

            doi = ""
            lid = str(rec.get("LID", "")).strip()
            if "[doi]" in lid.lower():
                doi = re.sub(r"\s*\[doi\]\s*", "", lid, flags=re.IGNORECASE).strip()
            if not doi:
                aids = rec.get("AID", [])
                if isinstance(aids, list):
                    for aid in aids:
                        aid_str = str(aid).strip()
                        if "[doi]" in aid_str.lower():
                            doi = re.sub(r"\s*\[doi\]\s*", "", aid_str, flags=re.IGNORECASE).strip()
                            break

            if title and abstract:
                papers.append({
                    "title": title,
                    "abstract": abstract,
                    "year": year,
                    "pmid": pmid,
                    "journal": journal,
                    "authors": authors,
                    "doi": doi,
                })

            if len(papers) >= max_results:
                break

        return papers
    except Exception:
        raise RuntimeError("No research data found")


def generate_research_paper(topic, papers):
    """
    Generate a structured literature review using Gemini with provided PubMed abstracts only.
    Uses a concise prompt and capped token output for faster response.
    """
    if not GEMINI_API_KEY:
        raise RuntimeError("GEMINI_API_KEY is not configured in .env")

    # Build a compact references block — only title + abstract (no full XML)
    references_block = "\n".join([
        f"{i+1}. {p['title']} ({p['year']}): {p['abstract'][:300]}"
        for i, p in enumerate(papers)
    ])

    # Short, focused prompt to reduce Gemini processing time
    prompt = (
        f"Write a concise literature review on '{topic}' using ONLY the studies below.\n"
        "Sections: Title, Abstract, Introduction, Methodology, Recent Findings, Discussion, Conclusion, References.\n"
        "Use only given titles in references. No invented facts. Be concise.\n\n"
        f"{references_block}"
    )

    payload = {
        "contents": [{"parts": [{"text": prompt}]}],
        "generationConfig": {
            "temperature": 0.2,
            "maxOutputTokens": 600,
        },
    }

    try:
        response = requests.post(
            f"https://generativelanguage.googleapis.com/v1beta/models/{GEMINI_MODEL}:generateContent?key={GEMINI_API_KEY}",
            headers={"Content-Type": "application/json"},
            json=payload,
            timeout=90,
        )
        if response.status_code == 429:
            raise RuntimeError("Gemini quota exceeded. Please check billing/usage and try again.")
        if response.status_code >= 400:
            raise RuntimeError("AI generation failed")

        data = response.json()
        candidates = data.get("candidates", [])
        if not candidates:
            raise RuntimeError("AI generation failed")

        parts = candidates[0].get("content", {}).get("parts", [])
        content = "\n".join([str(p.get("text", "")).strip() for p in parts if p.get("text")]).strip()
        if not content:
            raise RuntimeError("AI generation failed")

        return content
    except requests.exceptions.Timeout:
        raise RuntimeError("AI generation timed out. Please try again.")
    except requests.RequestException:
        raise RuntimeError("AI generation failed")


def extract_topic_from_uploaded_dna(filename, fasta_header):
    """Extract a probable biological topic from FASTA header or filename for NCBI lookup."""
    candidates = []
    if fasta_header:
        header = re.sub(r"^>+", "", fasta_header).strip()
        tokens = re.findall(r"[A-Za-z][A-Za-z0-9_-]{2,}", header)
        candidates.extend(tokens[:8])

    stem = os.path.splitext(os.path.basename(filename or ""))[0]
    file_tokens = re.findall(r"[A-Za-z][A-Za-z0-9_-]{2,}", stem)
    candidates.extend(file_tokens[:5])

    stop_words = {
        "fasta", "sequence", "sample", "dna", "gene", "upload", "untitled",
        "chromosome", "contig", "scaffold", "complete", "partial",
    }
    for token in candidates:
        t = token.strip("_- ")
        if not t:
            continue
        if t.lower() in stop_words:
            continue
        return t

    return stem or "DNA target"


def get_ncbi_context_for_topic(topic, max_results=3):
    """Fetch compact PubMed context from NCBI for simulation interpretability."""
    try:
        papers = search_pubmed_papers(topic, max_results=max_results)
    except Exception:
        papers = []

    return {
        "topic": topic,
        "count": len(papers),
        "papers": [
            {
                "title": p.get("title", ""),
                "year": p.get("year", "N/A"),
                "journal": p.get("journal", ""),
                "pmid": p.get("pmid", ""),
            }
            for p in papers
        ],
    }


def generate_ai_simulation_insight(sim_result, ncbi_context):
    """Generate short AI insight grounded in local simulation output plus NCBI references."""
    topic = ncbi_context.get("topic") or "DNA target"
    refs = ncbi_context.get("papers", [])[:3]
    bond = sim_result.get("bond_profile", {})

    if not GEMINI_API_KEY:
        best = sim_result.get("best_energy")
        avg = sim_result.get("average_energy")
        pose = sim_result.get("best_pose")
        protein_len = sim_result.get("protein_length")
        seq_len = sim_result.get("length")
        gc = sim_result.get("gc_content")
        return {
            "text": (
                f"For target topic {topic}, the uploaded DNA ({sim_result.get('protein_name')}) with length {seq_len} nt and GC content {gc}% "
                f"generated a best translated ORF protein length of {protein_len} aa, and ligand {sim_result.get('ligand_name')} produced a simulated docking profile with "
                f"best binding energy {best} kcal/mol at pose {pose} and average energy {avg} kcal/mol; the interaction pattern shows estimated hydrogen bonds "
                f"{bond.get('hydrogen_bonds')}, hydrophobic contacts {bond.get('hydrophobic_contacts')}, ionic interactions {bond.get('ionic_interactions')}, and "
                f"pi-stacking {bond.get('pi_stacking')} around predicted pocket {bond.get('predicted_binding_pocket')}, while NCBI evidence retrieval returned "
                f"{ncbi_context.get('count', 0)} supporting PubMed records for contextual interpretation."
            ),
            "status": "fallback",
            "error": "GEMINI_API_KEY is not configured.",
        }

    ref_lines = "\n".join([
        f"- {r.get('title', '')} ({r.get('year', 'N/A')}), PMID:{r.get('pmid', '')}"
        for r in refs
    ]) or "- No PubMed papers found for this topic"

    prompt = (
        "You are a bioinformatics assistant. Write exactly one detailed scientific paragraph (single paragraph, no bullet points, no headings, 120-180 words) "
        "strictly using the provided simulation metrics and NCBI references. Include DNA identity, ligand identity, core sequence properties, docking energies, "
        "bond-interaction profile, and NCBI evidence context. Do not invent wet-lab claims or extra data.\n\n"
        f"Topic: {topic}\n"
        f"DNA file: {sim_result.get('protein_name')}\n"
        f"Ligand file: {sim_result.get('ligand_name')}\n"
        f"Protein length (best ORF aa): {sim_result.get('protein_length')}\n"
        f"Best frame: {sim_result.get('best_frame')}\n"
        f"GC%: {sim_result.get('gc_content')}\n"
        f"Length(nt): {sim_result.get('length')}\n"
        f"Best energy: {sim_result.get('best_energy')}\n"
        f"Average energy: {sim_result.get('average_energy')}\n"
        f"Best pose: {sim_result.get('best_pose')}\n"
        f"Hydrogen bonds (estimated): {bond.get('hydrogen_bonds')}\n"
        f"Hydrophobic contacts (estimated): {bond.get('hydrophobic_contacts')}\n"
        f"Ionic interactions (estimated): {bond.get('ionic_interactions')}\n"
        f"Pi-stacking (estimated): {bond.get('pi_stacking')}\n"
        f"Predicted binding pocket: {bond.get('predicted_binding_pocket')}\n"
        "NCBI references:\n"
        f"{ref_lines}"
    )

    payload = {
        "contents": [{"parts": [{"text": prompt}]}],
        "generationConfig": {
            "temperature": 0.15,
            "maxOutputTokens": 320,
        },
    }

    last_error = "Unknown AI error"
    for api_version in ["v1beta", "v1"]:
        for model_name in _gemini_model_candidates():
            for _ in range(2):
                try:
                    response = requests.post(
                        f"https://generativelanguage.googleapis.com/{api_version}/models/{model_name}:generateContent?key={GEMINI_API_KEY}",
                        headers={"Content-Type": "application/json"},
                        json=payload,
                        timeout=45,
                    )
                    if response.status_code >= 400:
                        error_text = ""
                        try:
                            err_json = response.json()
                            error_text = str(err_json.get("error", {}).get("message", "")).strip()
                        except Exception:
                            error_text = response.text[:220]
                        last_error = (
                            f"Gemini {api_version} model {model_name} error {response.status_code}: "
                            f"{error_text or 'No error message'}"
                        )
                        continue

                    data = response.json()
                    candidates = data.get("candidates", [])
                    if not candidates:
                        last_error = f"Gemini {api_version} model {model_name} returned no candidates"
                        continue
                    parts = candidates[0].get("content", {}).get("parts", [])
                    text = "\n".join([str(p.get("text", "")).strip() for p in parts if p.get("text")]).strip()
                    if text:
                        one_paragraph = re.sub(r"\s*\n+\s*", " ", text).strip()
                        one_paragraph = re.sub(r"\s{2,}", " ", one_paragraph)
                        if len(one_paragraph) < 220:
                            one_paragraph = (
                                f"For topic {topic}, DNA file {sim_result.get('protein_name')} (length {sim_result.get('length')} nt, GC {sim_result.get('gc_content')}%) "
                                f"translated into a best ORF protein length of {sim_result.get('protein_length')} aa (frame {sim_result.get('best_frame')}), and ligand "
                                f"{sim_result.get('ligand_name')} produced a docking landscape with best binding energy {sim_result.get('best_energy')} kcal/mol at pose "
                                f"{sim_result.get('best_pose')} and average energy {sim_result.get('average_energy')} kcal/mol; the interaction profile estimates hydrogen bonds "
                                f"{bond.get('hydrogen_bonds')}, hydrophobic contacts {bond.get('hydrophobic_contacts')}, ionic interactions {bond.get('ionic_interactions')}, "
                                f"and pi-stacking {bond.get('pi_stacking')} near predicted pocket {bond.get('predicted_binding_pocket')}, while NCBI retrieval contributes "
                                f"{ncbi_context.get('count', 0)} related PubMed references for biological context."
                            )
                        return {
                            "text": one_paragraph,
                            "status": "ok",
                            "error": "",
                        }
                    last_error = f"Gemini {api_version} model {model_name} returned empty content"
                except requests.exceptions.Timeout:
                    last_error = f"Gemini {api_version} model {model_name} request timed out"
                except Exception as e:
                    last_error = f"Gemini {api_version} model {model_name} exception: {str(e)}"

    return {
        "text": (
            f"For topic {topic}, detailed AI generation is temporarily unavailable; based on current simulation and NCBI context, DNA file "
            f"{sim_result.get('protein_name')} with GC {sim_result.get('gc_content')}% and length {sim_result.get('length')} nt, alongside ligand "
            f"{sim_result.get('ligand_name')}, yielded best binding energy {sim_result.get('best_energy')} kcal/mol (pose {sim_result.get('best_pose')}) "
            f"and average {sim_result.get('average_energy')} kcal/mol with estimated interactions including hydrogen bonds {bond.get('hydrogen_bonds')}, "
            f"hydrophobic contacts {bond.get('hydrophobic_contacts')}, ionic interactions {bond.get('ionic_interactions')}, and pi-stacking {bond.get('pi_stacking')} "
            f"near pocket {bond.get('predicted_binding_pocket')}, while NCBI returned {ncbi_context.get('count', 0)} related PubMed records."
        ),
        "status": "unavailable",
        "error": last_error,
    }


def get_simulation_enrichment(topic, sim_result, max_results=2):
    """Cached enrichment so repeated runs do not block on NCBI/AI each time."""
    best_energy = sim_result.get("best_energy")
    avg_energy = sim_result.get("average_energy")
    cache_key = f"{str(topic).lower().strip()}|{best_energy}|{avg_energy}"
    now = time.time()

    cached = simulation_enrichment_cache.get(cache_key)
    if cached and (now - cached.get("timestamp", 0)) <= SIM_ENRICH_CACHE_TTL_SECONDS:
        return cached["data"]

    ncbi_context = get_ncbi_context_for_topic(topic, max_results=max_results)
    ai_result = generate_ai_simulation_insight(sim_result, ncbi_context)
    data = {
        "ncbi_context": ncbi_context,
        "ai_summary": ai_result.get("text", "AI insight unavailable."),
        "ai_status": ai_result.get("status", "unavailable"),
        "ai_error": ai_result.get("error", ""),
    }
    simulation_enrichment_cache[cache_key] = {
        "timestamp": now,
        "data": data,
    }
    return data


# ─────────────────────────────────────────────
# Flask Routes
# ─────────────────────────────────────────────

@app.route("/")
def index():
    """Serve the main single-page application."""
    return render_template("index.html")


@app.route("/api/search", methods=["POST"])
def api_search():
    """
    Accept a gene name and return top 3 NCBI results.
    POST /api/search  { "gene": "BRCA1" }
    """
    data = request.get_json()
    gene_name = data.get("gene", "").strip()

    if not gene_name:
        return jsonify({"error": "Please provide a gene name."}), 400

    # Basic input validation — allow only safe characters
    if not re.match(r"^[A-Za-z0-9\s\-_\.]+$", gene_name):
        return jsonify({"error": "Invalid gene name. Use only letters, numbers, spaces, or hyphens."}), 400

    try:
        results = search_gene(gene_name, max_results=3)
        if not results:
            return jsonify({"error": f"No results found for gene '{gene_name}'. Try a different name."}), 404
        return jsonify({"results": results})
    except RuntimeError as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/analyze", methods=["POST"])
def api_analyze():
    """
    Accept an NCBI ID and return full analysis results.
    POST /api/analyze  { "id": "NM_007294" }
    """
    data = request.get_json()
    ncbi_id = data.get("id", "").strip()

    if not ncbi_id:
        return jsonify({"error": "No sequence ID provided."}), 400

    try:
        record = fetch_sequence(ncbi_id)
        analysis = analyze_sequence(record)
        return jsonify(analysis)
    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    except RuntimeError as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/search-protein", methods=["POST"])
def api_search_protein():
    """
    Accept a protein name and return top 3 NCBI protein results.
    POST /api/search-protein  { "protein": "insulin" }
    """
    data = request.get_json()
    protein_name = data.get("protein", "").strip()

    if not protein_name:
        return jsonify({"error": "Please provide a protein name."}), 400

    if not re.match(r"^[A-Za-z0-9\s\-_\.]+$", protein_name):
        return jsonify({"error": "Invalid protein name. Use only letters, numbers, spaces, or hyphens."}), 400

    try:
        results = search_protein(protein_name, max_results=3)
        if not results:
            return jsonify({"error": f"No results found for protein '{protein_name}'. Try a different name."}), 404
        return jsonify({"results": results})
    except RuntimeError as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/analyze-protein", methods=["POST"])
def api_analyze_protein():
    """
    Accept an NCBI protein ID and return full protein analysis.
    POST /api/analyze-protein  { "id": "NP_000537" }
    """
    data = request.get_json()
    ncbi_id = data.get("id", "").strip()

    if not ncbi_id:
        return jsonify({"error": "No protein ID provided."}), 400

    try:
        record = fetch_protein_sequence(ncbi_id)
        analysis = analyze_protein_sequence(record)
        return jsonify(analysis)
    except RuntimeError as e:
        return jsonify({"error": str(e)}), 500


@app.route("/generate-paper", methods=["POST"])
@app.route("/api/research-paper", methods=["POST"])
def api_research_paper():
    """
    Generate a research paper using latest PubMed abstracts for a gene/protein topic.
    No AI/Gemini — builds the paper directly from PubMed data.
    POST /api/research-paper  { "topic": "BRCA1" }
    """
    data = request.get_json()
    topic = data.get("topic", "").strip() if data else ""
    analysis_context = data.get("analysis_context", {}) if data else {}
    if not isinstance(analysis_context, dict):
        analysis_context = {}

    if not topic:
        return jsonify({"error": "Please provide a gene or protein topic."}), 400

    if not re.match(r"^[A-Za-z0-9\s\-_\.\(\)]+$", topic):
        return jsonify({"error": "Invalid topic. Use only letters, numbers, spaces, hyphens, dots, or parentheses."}), 400

    # ── Check cache first ──
    topic_key = topic.lower().strip()
    context_key = json.dumps(analysis_context, sort_keys=True, default=str)
    cache_key = f"{topic_key}|{context_key}"
    cached = research_cache.get(cache_key)
    now = time.time()
    if cached and (now - cached.get("timestamp", 0)) <= RESEARCH_CACHE_TTL_SECONDS:
        return jsonify(cached["data"])

    # ── Fetch PubMed papers ──
    try:
        papers = search_pubmed_papers(topic, max_results=3)
    except Exception:
        papers = []

    if not papers:
        return jsonify({"error": "No research data found"}), 404

    # ── Build paper directly from PubMed abstracts (no Gemini) ──
    paper_text = build_paper_from_abstracts(topic, papers, analysis_context=analysis_context)

    response_payload = {
        "topic": topic,
        "paper": paper_text,
        "analysis_context_used": bool(analysis_context),
        "references": [
            {
                "title": p["title"],
                "year": p["year"],
                "journal": p.get("journal", ""),
                "pmid": p.get("pmid", ""),
                "doi": p.get("doi", ""),
            }
            for p in papers
        ],
    }
    research_cache[cache_key] = {
        "timestamp": now,
        "data": response_payload,
    }

    return jsonify(response_payload)


def _safe_float(value):
    try:
        return float(value)
    except Exception:
        return None


def _format_ieee_reference(idx, paper):
    authors = paper.get("authors", []) or []
    author_text = ""
    if authors:
        author_text = ", ".join(authors[:3]) + (", et al." if len(authors) > 3 else ".")
    else:
        author_text = "Unknown authors."

    title = str(paper.get("title", "Untitled")).strip().rstrip(".")
    journal = str(paper.get("journal", "")).strip()
    year = str(paper.get("year", "N/A")).strip()
    pmid = str(paper.get("pmid", "")).strip()
    doi = str(paper.get("doi", "")).strip()

    line = f"[{idx}] {author_text} \"{title},\""
    if journal:
        line += f" {journal},"
    if year and year != "N/A":
        line += f" {year}."
    else:
        line += " n.d."
    if doi:
        line += f" doi: {doi}."
    if pmid:
        line += f" PMID: {pmid}."
    return line


def _build_context_insights(analysis_context):
    dna = analysis_context.get("dna_analysis") if isinstance(analysis_context, dict) else None
    protein = analysis_context.get("protein_analysis") if isinstance(analysis_context, dict) else None
    drug = analysis_context.get("drug_simulation") if isinstance(analysis_context, dict) else None

    method_fragments = []
    results_fragments = []
    limitation_fragments = []

    if isinstance(dna, dict):
        seq_len = dna.get("sequence_length")
        gc = dna.get("gc_content")
        accession = dna.get("accession") or "N/A"
        if seq_len and gc is not None:
            method_fragments.append(
                f"The DNA analysis module processed accession {accession} with a nucleotide length of {seq_len} bp and computed a GC fraction of {gc}% using BioPython sequence utilities."
            )
            results_fragments.append(
                f"Sequence-level characterization identified a GC content of {gc}%, indicating target-dependent nucleotide composition relevant to coding stability and downstream translation."
            )

    if isinstance(protein, dict):
        prot_len = protein.get("sequence_length")
        mw = _safe_float(protein.get("molecular_weight"))
        p_i = protein.get("isoelectric_point")
        instability = _safe_float(protein.get("instability_index"))
        if prot_len:
            method_fragments.append(
                f"Protein physicochemical profiling was executed with ProtParam over a {prot_len} amino-acid sequence to derive molecular descriptors and stability proxies."
            )
        protein_result = ""
        if mw is not None:
            protein_result += f"Estimated molecular weight was {mw/1000:.2f} kDa"
        if p_i not in (None, "", "N/A"):
            protein_result += ("; " if protein_result else "") + f"theoretical pI was {p_i}"
        if instability is not None:
            protein_result += ("; " if protein_result else "") + f"instability index was {instability:.2f}"
        if protein_result:
            results_fragments.append(protein_result + ".")

    if isinstance(drug, dict):
        energy = _safe_float(drug.get("binding_energy"))
        disease = drug.get("disease") or "selected disease model"
        interaction = str(drug.get("interaction") or "interaction profile not reported").strip()
        demo_mode = bool(drug.get("demo_mode"))

        method_fragments.append(
            "Docking behavior was estimated through a simulated scoring layer that emulates ligand-target affinity trends without invoking external docking engines."
        )
        if energy is not None:
            results_fragments.append(
                f"In the {disease} model, the simulated lead binding energy was {energy:.2f} kcal/mol with interaction annotation indicating {interaction}."
            )
        if demo_mode:
            limitation_fragments.append(
                "Docking metrics were generated by a stochastic simulation module, therefore absolute energetic values should not be interpreted as experimentally calibrated free energies."
            )

    return method_fragments, results_fragments, limitation_fragments


def build_paper_from_abstracts(topic, papers, analysis_context=None):
    """Build a publication-style research paper from PubMed records and optional app analysis context."""
    analysis_context = analysis_context or {}
    years = [int(p["year"]) for p in papers if str(p.get("year", "")).isdigit()]
    year_window = f"{min(years)}-{max(years)}" if years else "recent years"

    condensed_findings = []
    for p in papers:
        abstract = str(p.get("abstract", "")).strip()
        first_sentence = re.split(r"(?<=[.!?])\s+", abstract)[0].strip() if abstract else ""
        if first_sentence:
            condensed_findings.append(first_sentence)

    method_ctx, result_ctx, limitation_ctx = _build_context_insights(analysis_context)

    title = (
        f"Integrative Web-Based Bioinformatics Pipeline for Sequence Characterization and "
        f"Simulated Molecular Docking in {topic}-Centered Therapeutic Discovery"
    )

    abstract = (
        f"{topic} remains a high-priority domain in translational bioinformatics, yet rapid integration of sequence analytics and candidate-level docking interpretation is often fragmented across independent tools. "
        f"This study presents a Flask-based web platform that unifies nucleotide profiling, protein property prediction, and simulated molecular docking within a single computational workflow. "
        f"The literature component was grounded in PubMed retrieval across {len(papers)} recent studies ({year_window}), followed by structured synthesis of study-level evidence. "
        f"The sequence module quantifies GC composition, codon usage, reverse complement, and translated peptide features, while the protein module estimates physicochemical descriptors, including molecular weight, theoretical isoelectric point, aromaticity, and instability index. "
        f"Docking output is represented as disease-aware simulated affinity scoring and interaction annotations to emulate early-stage screening behavior in a browser-native environment. "
        f"Application-linked outputs were incorporated into the evidence narrative to preserve consistency between literature context and computational findings. "
        f"The framework enables reproducible hypothesis refinement for target prioritization, while explicitly distinguishing simulated docking estimates from physics-based scoring engines."
    )

    introduction = (
        f"Computational drug discovery increasingly depends on integrated pipelines that couple sequence-level inference with ligand-target interaction assessment. "
        f"In many educational and exploratory environments, however, genomic analysis, protein interpretation, and docking-style outputs remain operationally isolated, which limits rapid iteration during hypothesis generation. "
        f"For {topic}, this fragmentation is particularly restrictive because target characterization and candidate screening are tightly interdependent stages. "
        f"Accordingly, the objective of this work was to design and evaluate a unified web-based system that performs DNA analytics, protein prediction, and simulated molecular docking, while maintaining direct traceability to current PubMed evidence."
    )

    methodology_core = (
        f"A structured PubMed query was executed against Title/Abstract fields using the term \"{topic}\", and the most recent {len(papers)} articles were selected for qualitative synthesis. "
        f"DNA analysis was implemented using BioPython sequence primitives to compute GC content, codon-frequency distributions, reverse-complement sequence generation, and codon-constrained translation. "
        f"Protein prediction logic applied ProtParam-derived physicochemical metrics, including molecular weight, theoretical pI, GRAVY, aromaticity, and secondary-structure fractions. "
        f"Drug discovery output was generated through simulated docking, wherein disease-conditioned random affinity scores and interaction summaries approximate early-stage virtual-screening readouts without invoking AutoDock-class engines. "
        f"System orchestration, API routing, and report generation were implemented in a Flask backend with client-side rendering for interactive exploration and export."
    )
    if method_ctx:
        methodology_core += " " + " ".join(method_ctx)

    results_text = (
        f"Evidence synthesis from the selected PubMed corpus demonstrated convergent emphasis on molecular mechanism, biomarker interpretation, and therapeutic modulation in {topic}. "
        + (" ".join(condensed_findings[:3]) if condensed_findings else "")
    )
    if result_ctx:
        results_text += " " + " ".join(result_ctx)
    else:
        results_text += (
            " The platform consistently produced sequence-derived descriptors and disease-specific simulated affinity values in ranges typically observed for exploratory in silico triage."
        )

    discussion = (
        f"The integrated workflow demonstrates that evidence-linked bioinformatics reporting can be operationalized in a lightweight web architecture without sacrificing analytical coherence. "
        f"A major strength is the coupling of literature retrieval with quantitative sequence and protein descriptors, enabling immediate contextual interpretation of computational outputs. "
        f"Another strength is transparent segregation of deterministic sequence analytics and probabilistic docking simulation, which reduces methodological ambiguity for end users. "
        f"The principal limitation is that docking outcomes are simulated rather than generated from force-field optimization or conformational sampling, thereby constraining external validity for lead ranking. "
        f"Additional limitations include dependence on abstract-level PubMed summaries, restricted corpus size, and potential bias introduced by recency-prioritized study selection."
    )
    if limitation_ctx:
        discussion += " " + " ".join(limitation_ctx)

    future_scope = (
        "Future development should integrate physics-based docking engines such as AutoDock Vina with explicit receptor and ligand preprocessing, binding-pocket definition, and pose-level scoring diagnostics. "
        "Model quality can be improved through AI-assisted ranking calibrated on benchmark complexes, uncertainty-aware affinity prediction, and retrieval-augmented literature synthesis incorporating full-text evidence rather than abstract-only data. "
        "A subsequent translational milestone would combine omics-derived target prioritization with active-learning driven compound suggestion to support iterative design cycles."
    )

    conclusion = (
        f"This study delivers a scientifically structured, web-based bioinformatics framework for {topic} that unifies DNA analysis, protein prediction, and simulated docking interpretation in a single reproducible environment. "
        "By grounding generated narratives in current PubMed records and app-level computational outputs, the system improves interpretability for early-stage exploratory research. "
        "Although docking values remain simulation-derived, the platform provides a robust foundation for migration toward fully validated structure-based drug discovery workflows."
    )

    references_lines = ["## References"]
    for idx, paper in enumerate(papers, start=1):
        references_lines.append(_format_ieee_reference(idx, paper))

    sections = [
        f"# Title\n{title}",
        f"## Abstract\n{abstract}",
        f"## Introduction\n{introduction}",
        f"## Methodology\n{methodology_core}",
        f"## Results\n{results_text}",
        f"## Discussion\n{discussion}",
        f"## Future Scope\n{future_scope}",
        f"## Conclusion\n{conclusion}",
        "\n".join(references_lines),
    ]
    return "\n\n".join(sections)


@app.route("/api/download-pdf", methods=["POST"])
def api_download_pdf():
    """
    Generate a PDF report with the full analysis data.
    POST /api/download-pdf  { full analysis JSON payload }
    """
    try:
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak
        from reportlab.lib import colors
        from reportlab.lib.units import mm
        from reportlab.lib.enums import TA_LEFT, TA_CENTER
    except ImportError:
        return jsonify({"error": "PDF library not installed. Run: pip install reportlab"}), 500

    data = request.get_json()
    if not data:
        return jsonify({"error": "No data provided."}), 400

    buffer = io.BytesIO()
    doc = SimpleDocTemplate(
        buffer, pagesize=A4,
        leftMargin=20*mm, rightMargin=20*mm,
        topMargin=20*mm, bottomMargin=20*mm
    )

    styles = getSampleStyleSheet()

    # Custom styles
    title_style = ParagraphStyle(
        'CustomTitle', parent=styles['Title'],
        fontSize=22, textColor=colors.HexColor('#1c1917'),
        spaceAfter=8, alignment=TA_CENTER
    )
    subtitle_style = ParagraphStyle(
        'CustomSubtitle', parent=styles['Normal'],
        fontSize=10, textColor=colors.HexColor('#57534e'),
        spaceAfter=20, alignment=TA_CENTER
    )
    heading_style = ParagraphStyle(
        'SectionHead', parent=styles['Heading2'],
        fontSize=14, textColor=colors.HexColor('#6366f1'),
        spaceBefore=18, spaceAfter=8
    )
    mono_style = ParagraphStyle(
        'MonoBlock', parent=styles['Normal'],
        fontName='Courier', fontSize=7, leading=10,
        textColor=colors.HexColor('#1c1917'),
        spaceAfter=10, wordWrap='CJK'
    )
    info_style = ParagraphStyle(
        'InfoText', parent=styles['Normal'],
        fontSize=10, textColor=colors.HexColor('#1c1917'),
        spaceAfter=6
    )

    story = []

    # Check if this is a protein report
    is_protein = data.get('_type') == 'protein'

    if is_protein:
        # ══════════ PROTEIN REPORT ══════════
        story.append(Paragraph("🧬 Protein Analysis Report", title_style))

        mw = data.get('molecular_weight', 'N/A')
        if mw != 'N/A':
            mw_display = f"{float(mw)/1000:.2f} kDa" if float(mw) > 1000 else f"{mw} Da"
        else:
            mw_display = 'N/A'

        story.append(Paragraph(
            f"Accession: {data.get('accession', 'N/A')} &nbsp;|&nbsp; "
            f"Length: {data.get('sequence_length', 'N/A')} aa &nbsp;|&nbsp; "
            f"Mol. Weight: {mw_display}",
            subtitle_style
        ))
        story.append(Paragraph(data.get('description', ''), info_style))
        story.append(Spacer(1, 10))

        # Summary Table
        summary_data = [
            ['Property', 'Value'],
            ['Accession', data.get('accession', 'N/A')],
            ['Sequence Length', f"{data.get('sequence_length', 'N/A')} aa"],
            ['Molecular Weight', mw_display],
            ['Isoelectric Point (pI)', str(data.get('isoelectric_point', 'N/A'))],
            ['Instability Index', f"{data.get('instability_index', 'N/A')} ({'Stable' if data.get('is_stable') else 'Unstable'})"],
            ['GRAVY Score', str(data.get('gravy', 'N/A'))],
            ['Aromaticity', str(data.get('aromaticity', 'N/A'))],
            ['Charge at pH 7.0', str(data.get('charge_at_ph7', 'N/A'))],
        ]

        t = Table(summary_data, colWidths=[80*mm, 80*mm])
        t.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#8b5cf6')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 9),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#e0e0e0')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f5f3ff')]),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('TOPPADDING', (0, 0), (-1, -1), 6),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
            ('LEFTPADDING', (0, 0), (-1, -1), 10),
        ]))
        story.append(t)
        story.append(Spacer(1, 12))

        # Secondary Structure
        sec = data.get('secondary_structure', {})
        if sec:
            story.append(Paragraph("🏛️ Secondary Structure Fraction", heading_style))
            sec_data = [
                ['Structure', 'Percentage'],
                ['α-Helix', f"{sec.get('helix', 0)}%"],
                ['Turn', f"{sec.get('turn', 0)}%"],
                ['β-Sheet', f"{sec.get('sheet', 0)}%"],
            ]
            st = Table(sec_data, colWidths=[80*mm, 80*mm])
            st.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#6366f1')),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, -1), 9),
                ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#e0e0e0')),
                ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f5f5f4')]),
                ('TOPPADDING', (0, 0), (-1, -1), 6),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
                ('LEFTPADDING', (0, 0), (-1, -1), 10),
            ]))
            story.append(st)
            story.append(Spacer(1, 12))

        # Full Protein Sequence
        full_seq = data.get('sequence_full', data.get('sequence_preview', 'N/A'))
        story.append(Paragraph("🧪 Full Protein Sequence", heading_style))
        story.append(Paragraph(f"Total amino acids: {len(full_seq)}", info_style))
        formatted_seq = format_seq_for_pdf(full_seq, 60)
        story.append(Paragraph(formatted_seq, mono_style))

        # Amino Acid Frequency
        aa_freq = data.get('amino_acid_frequency', {})
        aa_pct = data.get('amino_acid_percent', {})
        if aa_freq:
            story.append(Paragraph("📈 Amino Acid Frequency (Top 20)", heading_style))
            aa_data = [['Amino Acid', 'Count', 'Percentage']]
            for aa, count in aa_freq.items():
                pct = f"{aa_pct.get(aa, 0):.1f}%" if aa_pct else "—"
                aa_data.append([aa, str(count), pct])
            at = Table(aa_data, colWidths=[50*mm, 50*mm, 50*mm])
            at.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#8b5cf6')),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTNAME', (0, 1), (0, -1), 'Courier'),
                ('FONTSIZE', (0, 0), (-1, -1), 9),
                ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#e0e0e0')),
                ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f5f3ff')]),
                ('ALIGN', (1, 0), (-1, -1), 'CENTER'),
                ('TOPPADDING', (0, 0), (-1, -1), 5),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 5),
                ('LEFTPADDING', (0, 0), (-1, -1), 10),
            ]))
            story.append(at)

        filename = f"Protein_Report_{data.get('accession', 'unknown')}.pdf"

    else:
        # ══════════ DNA REPORT ══════════
        story.append(Paragraph("🧬 DNA & Protein Bioinformatics Analyzer — Full Report", title_style))
        story.append(Paragraph(
            f"Accession: {data.get('accession', 'N/A')} &nbsp;|&nbsp; "
            f"Length: {data.get('sequence_length', 'N/A')} bp &nbsp;|&nbsp; "
            f"GC Content: {data.get('gc_content', 'N/A')}%",
            subtitle_style
        ))
        story.append(Paragraph(data.get('description', ''), info_style))
        story.append(Spacer(1, 10))

        # Summary Table
        summary_data = [
            ['Property', 'Value'],
            ['Accession', data.get('accession', 'N/A')],
            ['Sequence Length', f"{data.get('sequence_length', 'N/A')} bp"],
            ['GC Content', f"{data.get('gc_content', 'N/A')}%"],
        ]
        nc = data.get('nucleotide_counts', {})
        if nc:
            summary_data.append(['Adenine (A)', str(nc.get('A', 0))])
            summary_data.append(['Thymine (T)', str(nc.get('T', 0))])
            summary_data.append(['Guanine (G)', str(nc.get('G', 0))])
            summary_data.append(['Cytosine (C)', str(nc.get('C', 0))])

        t = Table(summary_data, colWidths=[80*mm, 80*mm])
        t.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#6366f1')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 9),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#e0e0e0')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f5f5f4')]),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('TOPPADDING', (0, 0), (-1, -1), 6),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
            ('LEFTPADDING', (0, 0), (-1, -1), 10),
        ]))
        story.append(t)
        story.append(Spacer(1, 12))

        # Full DNA Sequence
        full_seq = data.get('sequence_full', data.get('sequence_preview', 'N/A'))
        story.append(Paragraph("🔬 Full DNA Sequence", heading_style))
        story.append(Paragraph(f"Total bases: {len(full_seq)}", info_style))
        formatted_seq = format_seq_for_pdf(full_seq, 60)
        story.append(Paragraph(formatted_seq, mono_style))

        # Full Reverse Complement
        full_rc = data.get('reverse_complement_full', data.get('reverse_complement', 'N/A'))
        story.append(Paragraph("🔄 Full Reverse Complement", heading_style))
        story.append(Paragraph(f"Total bases: {len(full_rc)}", info_style))
        formatted_rc = format_seq_for_pdf(full_rc, 60)
        story.append(Paragraph(formatted_rc, mono_style))

        # Full Protein Translation
        full_protein = data.get('protein_full', data.get('protein', 'N/A'))
        story.append(Paragraph("🧪 Full Protein Translation", heading_style))
        story.append(Paragraph(f"Total amino acids: {len(full_protein)}", info_style))
        formatted_protein = format_seq_for_pdf(full_protein, 60)
        story.append(Paragraph(formatted_protein, mono_style))

        # Codon Frequency
        codon_freq = data.get('codon_frequency', {})
        if codon_freq:
            story.append(Paragraph("📈 Codon Frequency (Top 20)", heading_style))
            codon_data = [['Codon', 'Count']]
            for codon, count in codon_freq.items():
                codon_data.append([codon, str(count)])
            ct = Table(codon_data, colWidths=[40*mm, 40*mm])
            ct.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#8b5cf6')),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTNAME', (0, 1), (0, -1), 'Courier'),
                ('FONTSIZE', (0, 0), (-1, -1), 9),
                ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#e0e0e0')),
                ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f5f3ff')]),
                ('ALIGN', (1, 0), (1, -1), 'CENTER'),
                ('TOPPADDING', (0, 0), (-1, -1), 5),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 5),
                ('LEFTPADDING', (0, 0), (-1, -1), 10),
            ]))
            story.append(ct)

        filename = f"DNA_Report_{data.get('accession', 'unknown')}.pdf"

    # Build PDF
    doc.build(story)
    buffer.seek(0)

    return send_file(
        buffer,
        mimetype='application/pdf',
        as_attachment=True,
        download_name=filename
    )


@app.route("/api/download-research-pdf", methods=["POST"])
def api_download_research_pdf():
    """
    Generate a PDF for AI-generated research paper output.
    POST /api/download-research-pdf  { topic, paper, references }
    """
    try:
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
        from reportlab.lib import colors
        from reportlab.lib.units import mm
        from reportlab.lib.enums import TA_LEFT, TA_CENTER
    except ImportError:
        return jsonify({"error": "PDF library not installed. Run: pip install reportlab"}), 500

    data = request.get_json()
    if not data:
        return jsonify({"error": "No data provided."}), 400

    topic = str(data.get("topic", "Research Topic")).strip() or "Research Topic"
    paper = str(data.get("paper", "")).strip()
    references = data.get("references", [])

    if not paper:
        return jsonify({"error": "No research paper content provided."}), 400

    buffer = io.BytesIO()
    doc = SimpleDocTemplate(
        buffer,
        pagesize=A4,
        leftMargin=20 * mm,
        rightMargin=20 * mm,
        topMargin=20 * mm,
        bottomMargin=20 * mm,
    )

    styles = getSampleStyleSheet()
    title_style = ParagraphStyle(
        "ResearchTitle",
        parent=styles["Title"],
        fontSize=20,
        textColor=colors.HexColor("#1c1917"),
        alignment=TA_CENTER,
        spaceAfter=8,
    )
    sub_style = ParagraphStyle(
        "ResearchSub",
        parent=styles["Normal"],
        fontSize=10,
        textColor=colors.HexColor("#57534e"),
        alignment=TA_CENTER,
        spaceAfter=16,
    )
    heading_style = ParagraphStyle(
        "ResearchHeading",
        parent=styles["Heading2"],
        fontSize=13,
        textColor=colors.HexColor("#6366f1"),
        alignment=TA_LEFT,
        spaceBefore=10,
        spaceAfter=6,
    )
    body_style = ParagraphStyle(
        "ResearchBody",
        parent=styles["Normal"],
        fontSize=10,
        textColor=colors.HexColor("#1c1917"),
        leading=14,
        alignment=TA_LEFT,
        spaceAfter=8,
    )

    story = [
        Paragraph("AI Research Paper", title_style),
        Paragraph(f"Topic: {topic}", sub_style),
    ]

    lines = [line.strip() for line in paper.splitlines() if line.strip()]
    headings = {
        "title",
        "abstract",
        "introduction",
        "methodology",
        "results",
        "recent findings",
        "discussion",
        "future scope",
        "conclusion",
        "references",
    }

    for line in lines:
        normalized = re.sub(r"^[#\-*\s]+", "", line).strip()
        normalized_key = normalized.rstrip(":").strip().lower()

        if normalized_key in headings:
            story.append(Paragraph(normalized.rstrip(":"), heading_style))
        elif any(normalized.lower().startswith(h + ":") for h in headings):
            label, value = normalized.split(":", 1)
            story.append(Paragraph(label.strip(), heading_style))
            if value.strip():
                story.append(Paragraph(value.strip().replace("\n", "<br/>"), body_style))
        else:
            story.append(Paragraph(normalized.replace("\n", "<br/>"), body_style))

    if references:
        story.append(Paragraph("References", heading_style))
        for idx, ref in enumerate(references, start=1):
            title = str(ref.get("title", "")).strip()
            year = str(ref.get("year", "N/A")).strip()
            if title:
                story.append(Paragraph(f"{idx}. {title} ({year})", body_style))

    doc.build(story)
    buffer.seek(0)

    safe_topic = re.sub(r"[^A-Za-z0-9_\-]+", "_", topic)[:40] or "research"
    filename = f"Research_Paper_{safe_topic}.pdf"

    return send_file(
        buffer,
        mimetype="application/pdf",
        as_attachment=True,
        download_name=filename,
    )


def format_seq_for_pdf(seq, line_len=60):
    """
    Format a sequence string into fixed-width lines for PDF display.
    Adds line numbers and spaces every 10 chars for readability.
    """
    if not seq:
        return "(No data available)"
    lines = []
    for i in range(0, len(seq), line_len):
        chunk = seq[i:i+line_len]
        # Add spaces every 10 chars
        spaced = ' '.join([chunk[j:j+10] for j in range(0, len(chunk), 10)])
        line_num = str(i+1).rjust(6)
        lines.append(f"{line_num}  {spaced}")
    return '<br/>'.join(lines)


# ─────────────────────────────────────────────
# Drug Discovery Simulation Routes
# ─────────────────────────────────────────────

@app.route("/drug-discovery")
def drug_discovery_ui():
    """Serve the standalone Drug Discovery Simulation page."""
    return render_template("drug_simulation.html")

@app.route("/api/drug-simulation", methods=["POST"])
def api_drug_simulation():
    """
    Handle DNA-driven docking simulation request.
    Accepts multipart/form-data:
      dna_file: required (.fasta, .fa, .fna, .txt, .seq)
      ligand_file: optional .sdf or .mol
      demo_project: optional built-in demo key for no-upload run
      enable_enrichment: optional checkbox for NCBI + AI insights (slower)
    """
    try:
        dna_file = request.files.get("dna_file") or request.files.get("protein_file")
        ligand_file = request.files.get("ligand_file")
        demo_project = (request.form.get("demo_project") or "").strip()

        has_protein = bool(dna_file and dna_file.filename)
        has_ligand = bool(ligand_file and ligand_file.filename)
        enable_enrichment = (request.form.get("enable_enrichment") or "").lower() in {"on", "true", "1", "yes"}

        protein_name = dna_file.filename if has_protein else None
        ligand_name = ligand_file.filename if has_ligand else None

        use_demo = bool((not has_protein) and demo_project)
        demo_cfg = DEMO_PROJECTS.get(demo_project) if use_demo else None

        if use_demo and not demo_cfg:
            return jsonify({"error": "Invalid demo project selected."}), 400

        if not has_protein and not use_demo:
            return jsonify({"error": "DNA file is required, or choose a demo project."}), 400

        allowed_dna = {".fasta", ".fa", ".fna", ".txt", ".seq"}
        if has_protein and not any(protein_name.lower().endswith(ext) for ext in allowed_dna):
            return jsonify({"error": "DNA file format not supported. Use .fasta, .fa, .fna, .txt, or .seq"}), 400
        if has_ligand and not (ligand_name.lower().endswith(".sdf") or ligand_name.lower().endswith(".mol")):
            return jsonify({"error": "Ligand file must be .sdf or .mol"}), 400

        protein_sequence = ""
        fasta_header = ""
        if use_demo and demo_cfg:
            protein_name = demo_cfg["dna_name"]
            ligand_name = ligand_name or demo_cfg["ligand_name"]
            fasta_header = demo_cfg.get("fasta_header", "")
            protein_sequence = str(demo_cfg.get("dna_sequence", "")).upper()
            has_protein = True
            has_ligand = bool(ligand_name)
        elif has_protein:
            raw = dna_file.read().decode("utf-8", errors="ignore")
            header_line = next((line.strip() for line in raw.splitlines() if line.strip().startswith(">")), "")
            fasta_header = header_line
            cleaned_lines = [line.strip() for line in raw.splitlines() if not line.strip().startswith(">")]
            cleaned_text = "".join(cleaned_lines).upper()
            protein_sequence = "".join(ch for ch in cleaned_text if ch in {"A", "C", "G", "T", "N"})

            if not protein_sequence:
                return jsonify({"error": "Empty DNA sequence. File does not contain A/C/G/T bases."}), 400

        result = simulate_docking(
            has_protein=has_protein,
            has_ligand=has_ligand,
            protein_name=protein_name,
            ligand_name=ligand_name,
            protein_sequence=protein_sequence,
        )

        topic = extract_topic_from_uploaded_dna(protein_name, fasta_header)

        result["dna_name"] = protein_name
        result["ligand_name"] = ligand_name or "No ligand uploaded"
        result["demo_project"] = demo_project if use_demo else ""
        result["used_demo"] = use_demo
        result["extracted_topic"] = topic
        result["enrichment_enabled"] = enable_enrichment

        if enable_enrichment:
            enrichment = get_simulation_enrichment(topic, result, max_results=2)
            result["ncbi_context"] = enrichment.get("ncbi_context", {})
            result["ai_summary"] = enrichment.get("ai_summary", "AI insight unavailable.")
            result["ai_status"] = enrichment.get("ai_status", "unavailable")
            result["ai_error"] = enrichment.get("ai_error", "")
        else:
            result["ncbi_context"] = {"topic": topic, "count": 0, "papers": []}
            result["ai_summary"] = "Fast mode enabled: NCBI + AI enrichment skipped for quicker loading."
            result["ai_status"] = "skipped"
            result["ai_error"] = ""

        return jsonify(result), 200
    except ValueError as e:
        return jsonify({"error": str(e), "message": "Invalid biological input"}), 400
    except Exception as e:
        return jsonify({"error": str(e), "message": "Failed to run drug simulation"}), 500

if __name__ == "__main__":
    # Run in debug mode on port 5000 for local development
    # threaded=True allows concurrent requests so the UI doesn't freeze
    app.run(debug=True, port=5000, threaded=True)
