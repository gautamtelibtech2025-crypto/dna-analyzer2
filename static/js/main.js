/**
 * main.js – DNA & Protein Bioinformatics Analyzer
 * Handles gene/protein search, result selection, sequence analysis, and rendering.
 * Uses Chart.js for charts.
 */

// ─── State ────────────────────────────────────────────────────
let gcChart = null;
let aaChart = null;
let secChart = null;
let lastAnalysisData = null;
let lastProteinData = null;
let lastResearchPaperData = null;
let currentMode = 'dna'; // 'dna' or 'protein'

// ─── DOM Helpers ──────────────────────────────────────────────
const $ = (id) => document.getElementById(id);

function show(id) { $(id).classList.remove("hidden"); }
function hide(id) { $(id).classList.add("hidden"); }
function setText(id, val) { $(id).textContent = val; }

// ─── Mode Toggle ─────────────────────────────────────────────
function setMode(mode) {
  currentMode = mode;

  // Toggle button active states
  $("mode-dna").classList.toggle("active", mode === 'dna');
  $("mode-protein").classList.toggle("active", mode === 'protein');

  // Update placeholder and quick tags
  if (mode === 'dna') {
    $("gene-input").placeholder = "e.g. BRCA1, TP53, EGFR, MYC";
    show("dna-quick-tags");
    hide("protein-quick-tags");
  } else {
    $("gene-input").placeholder = "e.g. insulin, hemoglobin, p53, collagen";
    hide("dna-quick-tags");
    show("protein-quick-tags");
  }

  // Hide any existing results
  hide("error-msg");
  hide("results-selection");
  hide("analysis-section");
  hide("protein-analysis-section");
  hide("fetch-spinner");
}

// ─── Quick Search (from tag buttons) ─────────────────────────
function quickSearch(gene) {
  $("gene-input").value = gene;
  searchGene();
}

// ─── Entry Point: User clicks "Analyze" or presses Enter ────
async function searchGene() {
  const input = $("gene-input").value.trim();
  if (!input) {
    $("gene-input").focus();
    return;
  }

  // Reset UI
  hide("error-msg");
  hide("results-selection");
  hide("analysis-section");
  hide("protein-analysis-section");
  hide("fetch-spinner");
  setButtonLoading(true);

  try {
    let data;
    if (currentMode === 'dna') {
      data = await postJSON("/api/search", { gene: input });
    } else {
      data = await postJSON("/api/search-protein", { protein: input });
    }

    if (data.error) {
      showError(data.error);
      return;
    }

    renderResultSelection(data.results);
  } catch (err) {
    showError("Network error – is the Flask server running?");
  } finally {
    setButtonLoading(false);
  }
}

// ─── Render the 3 search result cards ─────────────────────────
function renderResultSelection(results) {
  const container = $("result-cards");
  container.innerHTML = "";

  const unit = currentMode === 'dna' ? 'bp' : 'aa';

  results.forEach((r, idx) => {
    const div = document.createElement("div");
    div.className = "result-item";
    div.setAttribute("role", "button");
    div.setAttribute("tabindex", "0");
    div.style.animationDelay = `${idx * 0.08}s`;
    div.innerHTML = `
      <div class="result-num">${idx + 1}</div>
      <div class="result-info">
        <div class="result-title" title="${escHtml(r.title)}">${escHtml(r.title)}</div>
        <div class="result-meta">ID: ${r.id} &nbsp;&middot;&nbsp; Length: ${Number(r.length).toLocaleString()} ${unit}</div>
      </div>
      <div class="result-arrow">›</div>
    `;

    if (currentMode === 'dna') {
      div.addEventListener("click", () => analyzeSequence(r.id, r.title));
      div.addEventListener("keydown", (e) => {
        if (e.key === "Enter") analyzeSequence(r.id, r.title);
      });
    } else {
      div.addEventListener("click", () => analyzeProtein(r.id, r.title));
      div.addEventListener("keydown", (e) => {
        if (e.key === "Enter") analyzeProtein(r.id, r.title);
      });
    }

    container.appendChild(div);
  });

  show("results-selection");
}

// ─── Fetch full DNA analysis for the selected sequence ────────
async function analyzeSequence(ncbiId, title) {
  hide("results-selection");
  hide("error-msg");
  show("fetch-spinner");

  try {
    const data = await postJSON("/api/analyze", { id: ncbiId, title });

    if (data.error) {
      hide("fetch-spinner");
      showError(data.error);
      show("results-selection");
      return;
    }

    renderAnalysis(data);
  } catch (err) {
    showError("Network error during sequence fetch.");
    show("results-selection");
  } finally {
    hide("fetch-spinner");
  }
}

// ─── Fetch full Protein analysis ──────────────────────────────
async function analyzeProtein(ncbiId, title) {
  hide("results-selection");
  hide("error-msg");
  show("fetch-spinner");

  try {
    const data = await postJSON("/api/analyze-protein", { id: ncbiId, title });

    if (data.error) {
      hide("fetch-spinner");
      showError(data.error);
      show("results-selection");
      return;
    }

    renderProteinAnalysis(data);
  } catch (err) {
    showError("Network error during protein fetch.");
    show("results-selection");
  } finally {
    hide("fetch-spinner");
  }
}

// ─── Populate the DNA analysis results section ────────────────
function renderAnalysis(d) {
  lastAnalysisData = d;

  setText("res-accession", d.accession);
  setText("res-length", Number(d.sequence_length).toLocaleString() + " bp");
  setText("res-gc", d.gc_content + "%");
  setText("res-description", d.description);

  setText("res-sequence", formatDNA(d.sequence_preview));
  if (d.sequence_full) setText("res-sequence-full", formatDNA(d.sequence_full));

  setText("res-revcomp", formatDNA(d.reverse_complement));
  if (d.reverse_complement_full) setText("res-revcomp-full", formatDNA(d.reverse_complement_full));

  setText("res-protein", d.protein || "(No protein translated – check reading frame)");
  if (d.protein_full) setText("res-protein-full", d.protein_full || "(No protein data)");

  renderNucleotideChart(d.gc_content, d.nucleotide_counts);

  const at = (100 - d.gc_content).toFixed(2);
  $("gc-summary").textContent =
    `GC: ${d.gc_content}%  ·  AT: ${at}%  (${Number(d.sequence_length).toLocaleString()} total bases)`;

  renderCodonTable(d.codon_frequency);

  hide("full-sequence-wrap");
  hide("full-revcomp-wrap");
  hide("full-protein-wrap");

  show("analysis-section");
  hide("protein-analysis-section");

  setTimeout(() => $("analysis-section").scrollIntoView({ behavior: "smooth", block: "start" }), 80);
}

// ─── Populate the Protein analysis results section ────────────
function renderProteinAnalysis(d) {
  lastProteinData = d;

  // Banner
  setText("prot-accession", d.accession);
  setText("prot-length", Number(d.sequence_length).toLocaleString() + " aa");

  // Molecular weight formatting
  if (d.molecular_weight !== "N/A") {
    const mw = Number(d.molecular_weight);
    if (mw > 1000) {
      setText("prot-mw", (mw / 1000).toFixed(2) + " kDa");
    } else {
      setText("prot-mw", mw.toFixed(2) + " Da");
    }
  } else {
    setText("prot-mw", "N/A");
  }

  setText("prot-description", d.description);

  // Protein sequence preview
  setText("prot-sequence", formatProtein(d.sequence_preview));
  if (d.sequence_full) setText("prot-sequence-full", formatProtein(d.sequence_full));

  // Physicochemical properties
  setText("prot-pi", d.isoelectric_point);
  setText("prot-instability", d.instability_index);
  setText("prot-gravy", d.gravy);
  setText("prot-aromaticity", d.aromaticity);
  setText("prot-charge", d.charge_at_ph7);

  // Stability tag
  const stabTag = $("prot-stability-tag");
  if (d.is_stable === true) {
    stabTag.textContent = "STABLE";
    stabTag.className = "prop-tag tag-stable";
  } else if (d.is_stable === false) {
    stabTag.textContent = "UNSTABLE";
    stabTag.className = "prop-tag tag-unstable";
  } else {
    stabTag.textContent = "";
  }

  // Amino acid composition chart
  renderAAChart(d.amino_acid_count);

  // Secondary structure chart
  renderSecStructureChart(d.secondary_structure);

  // Secondary structure summary
  if (d.secondary_structure) {
    const ss = d.secondary_structure;
    $("sec-summary").textContent =
      `Helix: ${ss.helix}%  ·  Turn: ${ss.turn}%  ·  Sheet: ${ss.sheet}%`;
  }

  // Amino acid frequency table
  renderAATable(d.amino_acid_frequency, d.amino_acid_percent);

  // Reset full views
  hide("full-protseq-wrap");

  hide("analysis-section");
  show("protein-analysis-section");

  setTimeout(() => $("protein-analysis-section").scrollIntoView({ behavior: "smooth", block: "start" }), 80);
}

// ─── Amino Acid Composition Chart (Bar) ───────────────────────
function renderAAChart(aaCounts) {
  const ctx = $("aa-chart").getContext("2d");

  if (aaChart) aaChart.destroy();

  if (!aaCounts || Object.keys(aaCounts).length === 0) return;

  // Amino acid full names
  const aaNames = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'
  };

  const sorted = Object.entries(aaCounts).sort((a, b) => b[1] - a[1]);
  const labels = sorted.map(([aa]) => `${aa} (${aaNames[aa] || aa})`);
  const values = sorted.map(([, count]) => count);

  // Generate gradient colors
  const colors = sorted.map((_, i) => {
    const hue = (i * 18) + 200;
    return `hsla(${hue % 360}, 70%, 55%, 0.8)`;
  });

  aaChart = new Chart(ctx, {
    type: "bar",
    data: {
      labels,
      datasets: [{
        data: values,
        backgroundColor: colors,
        borderColor: colors.map(c => c.replace('0.8', '1')),
        borderWidth: 1,
        borderRadius: 4,
      }],
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      plugins: {
        legend: { display: false },
        tooltip: {
          callbacks: {
            label: (ctx) => ` Count: ${ctx.raw.toLocaleString()}`,
          },
        },
      },
      scales: {
        x: {
          grid: { display: false },
          ticks: { font: { size: 9, family: "JetBrains Mono" }, maxRotation: 45 },
        },
        y: {
          grid: { color: "rgba(0,0,0,.04)" },
          ticks: { font: { size: 10, family: "Inter" } },
        },
      },
    },
  });
}

// ─── Secondary Structure Chart (Doughnut) ─────────────────────
function renderSecStructureChart(secStructure) {
  const ctx = $("sec-chart").getContext("2d");

  if (secChart) secChart.destroy();

  if (!secStructure) return;

  secChart = new Chart(ctx, {
    type: "doughnut",
    data: {
      labels: ["α-Helix", "Turn", "β-Sheet"],
      datasets: [{
        data: [secStructure.helix, secStructure.turn, secStructure.sheet],
        backgroundColor: [
          "rgba(99, 102, 241, .8)",   // Indigo for helix
          "rgba(245, 158, 11, .8)",   // Amber for turn
          "rgba(16, 185, 129, .8)",   // Emerald for sheet
        ],
        borderColor: "#fff",
        borderWidth: 2.5,
        hoverOffset: 6,
      }],
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      cutout: "62%",
      plugins: {
        legend: {
          position: "bottom",
          labels: {
            padding: 12,
            usePointStyle: true,
            pointStyleWidth: 8,
            font: { size: 11, family: "Inter" },
          },
        },
        tooltip: {
          callbacks: {
            label: (ctx) => ` ${ctx.label}: ${ctx.raw}%`,
          },
        },
      },
    },
  });
}

// ─── Amino Acid Frequency Table ───────────────────────────────
function renderAATable(freqObj, percentObj) {
  const tbody = $("aa-tbody");
  tbody.innerHTML = "";

  if (!freqObj || Object.keys(freqObj).length === 0) {
    tbody.innerHTML =
      '<tr><td colspan="4" style="color:var(--text-3);text-align:center;padding:20px">No amino acid data available</td></tr>';
    return;
  }

  const aaNames = {
    'A': 'Alanine', 'R': 'Arginine', 'N': 'Asparagine', 'D': 'Aspartate',
    'C': 'Cysteine', 'E': 'Glutamate', 'Q': 'Glutamine', 'G': 'Glycine',
    'H': 'Histidine', 'I': 'Isoleucine', 'L': 'Leucine', 'K': 'Lysine',
    'M': 'Methionine', 'F': 'Phenylalanine', 'P': 'Proline', 'S': 'Serine',
    'T': 'Threonine', 'W': 'Tryptophan', 'Y': 'Tyrosine', 'V': 'Valine'
  };

  const entries = Object.entries(freqObj);
  const maxCount = entries[0][1];

  entries.forEach(([aa, count]) => {
    const pct = ((count / maxCount) * 100).toFixed(0);
    const percent = percentObj && percentObj[aa] ? percentObj[aa].toFixed(1) + "%" : "—";
    const tr = document.createElement("tr");
    tr.innerHTML = `
      <td class="codon-td-codon">${aa} <span style="color:var(--text-3);font-family:Inter;font-weight:400;font-size:.72rem">${aaNames[aa] || ''}</span></td>
      <td>${count.toLocaleString()}</td>
      <td>${percent}</td>
      <td class="codon-bar-cell">
        <div class="codon-bar-wrap">
          <div class="codon-bar-fill protein-bar-fill" style="width:${pct}%"></div>
        </div>
      </td>
    `;
    tbody.appendChild(tr);
  });
}

// ─── Nucleotide Composition Chart (Chart.js) ──────────────────
function renderNucleotideChart(gcPercent, nucleotideCounts) {
  const ctx = $("gc-chart").getContext("2d");
  if (gcChart) gcChart.destroy();

  if (nucleotideCounts) {
    const { A, T, G, C } = nucleotideCounts;
    gcChart = new Chart(ctx, {
      type: "doughnut",
      data: {
        labels: ["Adenine (A)", "Thymine (T)", "Guanine (G)", "Cytosine (C)"],
        datasets: [{
          data: [A, T, G, C],
          backgroundColor: [
            "rgba(59, 130, 246, .8)",
            "rgba(245, 158, 11, .8)",
            "rgba(16, 185, 129, .8)",
            "rgba(139, 92, 246, .8)",
          ],
          borderColor: "#fff",
          borderWidth: 2.5,
          hoverOffset: 6,
        }],
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        cutout: "62%",
        plugins: {
          legend: {
            position: "bottom",
            labels: { padding: 12, usePointStyle: true, pointStyleWidth: 8, font: { size: 11, family: "Inter" } },
          },
          tooltip: {
            callbacks: {
              label: (ctx) => {
                const total = ctx.dataset.data.reduce((a, b) => a + b, 0);
                const pct = ((ctx.raw / total) * 100).toFixed(1);
                return ` ${ctx.label}: ${ctx.raw.toLocaleString()} (${pct}%)`;
              },
            },
          },
        },
      },
    });
  } else {
    const atPercent = parseFloat((100 - gcPercent).toFixed(2));
    gcChart = new Chart(ctx, {
      type: "bar",
      data: {
        labels: ["GC Content", "AT Content"],
        datasets: [{
          data: [gcPercent, atPercent],
          backgroundColor: ["rgba(16,185,129,.8)", "rgba(59,130,246,.7)"],
          borderColor: ["rgba(16,185,129,1)", "rgba(59,130,246,1)"],
          borderWidth: 1.5,
          borderRadius: 8,
        }],
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
          legend: { display: false },
          tooltip: { callbacks: { label: (ctx) => ` ${ctx.raw}%` } },
        },
        scales: {
          x: { grid: { display: false }, ticks: { font: { size: 11, family: "Inter" } } },
          y: {
            min: 0, max: 100,
            grid: { color: "rgba(0,0,0,.04)" },
            ticks: { callback: (v) => v + "%", font: { size: 10, family: "Inter" }, stepSize: 25 },
          },
        },
      },
    });
  }
}

// ─── Codon Frequency Table ────────────────────────────────────
function renderCodonTable(freqObj) {
  const tbody = $("codon-tbody");
  tbody.innerHTML = "";

  if (!freqObj || Object.keys(freqObj).length === 0) {
    tbody.innerHTML =
      '<tr><td colspan="3" style="color:var(--text-3);text-align:center;padding:20px">No codon data available</td></tr>';
    return;
  }

  const entries = Object.entries(freqObj);
  const maxCount = entries[0][1];

  entries.forEach(([codon, count]) => {
    const pct = ((count / maxCount) * 100).toFixed(0);
    const tr = document.createElement("tr");
    tr.innerHTML = `
      <td class="codon-td-codon">${codon}</td>
      <td>${count.toLocaleString()}</td>
      <td class="codon-bar-cell">
        <div class="codon-bar-wrap">
          <div class="codon-bar-fill" style="width:${pct}%"></div>
        </div>
      </td>
    `;
    tbody.appendChild(tr);
  });
}

// ─── Utility Functions ────────────────────────────────────────

/** Format DNA into 10-base blocks separated by spaces, 6 blocks per line */
function formatDNA(seq) {
  if (!seq) return "";
  const blocks = seq.match(/.{1,10}/g) || [];
  const lines = [];
  for (let i = 0; i < blocks.length; i += 6) {
    lines.push(blocks.slice(i, i + 6).join(" "));
  }
  return lines.join("\n");
}

/** Format protein into 10-aa blocks */
function formatProtein(seq) {
  if (!seq) return "";
  const blocks = seq.match(/.{1,10}/g) || [];
  const lines = [];
  for (let i = 0; i < blocks.length; i += 6) {
    lines.push(blocks.slice(i, i + 6).join(" "));
  }
  return lines.join("\n");
}

/** Escape HTML characters to prevent XSS */
function escHtml(str) {
  return str
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;");
}

/** POST JSON helper */
async function postJSON(url, body, timeoutMs = 120000) {
  const controller = new AbortController();
  const timer = setTimeout(() => controller.abort(), timeoutMs);

  try {
    const res = await fetch(url, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
      signal: controller.signal,
    });
    const json = await res.json();
    return json;
  } catch (err) {
    if (err.name === "AbortError") {
      return { error: "Request timeout. The server is taking too long. Please try again." };
    }
    throw err;
  } finally {
    clearTimeout(timer);
  }
}

/** Show the error message banner */
function showError(msg) {
  const el = $("error-msg");
  el.textContent = "⚠ " + msg;
  show("error-msg");
}

/** Toggle search button loading state */
function setButtonLoading(loading) {
  $("analyze-btn").disabled = loading;
  if (loading) {
    hide("btn-label");
    show("btn-spinner");
  } else {
    show("btn-label");
    hide("btn-spinner");
  }
}

/** Reset app */
function resetApp() {
  hide("results-selection");
  hide("analysis-section");
  hide("protein-analysis-section");
  hide("error-msg");
  hide("full-sequence-wrap");
  hide("full-revcomp-wrap");
  hide("full-protein-wrap");
  hide("full-protseq-wrap");
  lastAnalysisData = null;
  lastProteinData = null;
  lastResearchPaperData = null;
  hide("research-error");
  hide("research-output");
  if ($("research-input")) {
    $("research-input").value = "";
  }
  $("gene-input").value = "";
  $("gene-input").focus();
  window.scrollTo({ top: 0, behavior: "smooth" });
}

// ─── Toggle Full Sequence View ────────────────────────────────
function toggleFullView(type) {
  const map = {
    'sequence': ['full-sequence-wrap', 'toggle-seq-btn'],
    'revcomp': ['full-revcomp-wrap', 'toggle-revcomp-btn'],
    'protein': ['full-protein-wrap', 'toggle-protein-btn'],
    'protseq': ['full-protseq-wrap', 'toggle-protseq-btn'],
  };
  const [wrapId, btnId] = map[type] || map['sequence'];
  const wrap = $(wrapId);
  const btn = $(btnId);

  if (wrap.classList.contains('hidden')) {
    wrap.classList.remove('hidden');
    wrap.style.animation = 'fadeUp .35s var(--ease)';
    btn.querySelector('span').textContent = 'Hide Full';
    btn.classList.add('active');
  } else {
    wrap.classList.add('hidden');
    btn.querySelector('span').textContent = 'View Full';
    btn.classList.remove('active');
  }
}

// ─── Download DNA PDF Report ──────────────────────────────────
async function downloadPDF() {
  if (!lastAnalysisData) {
    showError("No analysis data available. Please analyze a gene first.");
    return;
  }

  const btn = $("download-pdf-btn");
  const label = $("pdf-btn-label");
  const spinner = $("pdf-btn-spinner");

  btn.disabled = true;
  label.classList.add('hidden');
  spinner.classList.remove('hidden');

  try {
    const res = await fetch("/api/download-pdf", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(lastAnalysisData),
    });

    if (!res.ok) {
      const err = await res.json();
      showError(err.error || "Failed to generate PDF.");
      return;
    }

    const blob = await res.blob();
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `DNA_Report_${lastAnalysisData.accession || 'unknown'}.pdf`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  } catch (err) {
    showError("Network error while generating PDF.");
  } finally {
    btn.disabled = false;
    label.classList.remove('hidden');
    spinner.classList.add('hidden');
  }
}

// ─── Download Protein PDF Report ──────────────────────────────
async function downloadProteinPDF() {
  if (!lastProteinData) {
    showError("No protein data available. Please analyze a protein first.");
    return;
  }

  const btn = $("download-protein-pdf-btn");
  const label = $("protein-pdf-btn-label");
  const spinner = $("protein-pdf-btn-spinner");

  btn.disabled = true;
  label.classList.add('hidden');
  spinner.classList.remove('hidden');

  try {
    const res = await fetch("/api/download-pdf", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ ...lastProteinData, _type: "protein" }),
    });

    if (!res.ok) {
      const err = await res.json();
      showError(err.error || "Failed to generate PDF.");
      return;
    }

    const blob = await res.blob();
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `Protein_Report_${lastProteinData.accession || 'unknown'}.pdf`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  } catch (err) {
    showError("Network error while generating PDF.");
  } finally {
    btn.disabled = false;
    label.classList.remove('hidden');
    spinner.classList.add('hidden');
  }
}

// ─── AI Research Paper Generator ─────────────────────────────
async function generateResearchPaper() {
  const topicInput = $("research-input");
  if (!topicInput) return;

  const topic = topicInput.value.trim();
  if (!topic) {
    topicInput.focus();
    return;
  }

  hide("research-error");
  hide("research-output");
  setResearchButtonLoading(true);

  // Show a progress message that updates while waiting
  const progressMessages = [
    "Searching PubMed for recent papers…",
    "Extracting evidence from PubMed abstracts…",
    "Integrating sequence and protein analytics…",
    "Generating structured scientific manuscript…",
  ];
  let msgIdx = 0;
  const btnLabel = $("research-btn-label");
  const progressTimer = setInterval(() => {
    msgIdx++;
    if (msgIdx < progressMessages.length && btnLabel) {
      btnLabel.textContent = progressMessages[msgIdx];
      btnLabel.classList.remove("hidden");
    }
  }, 8000);

  try {
    const analysisContext = buildResearchContext();

    // 2 min timeout — PubMed + Gemini + buffer
    const data = await postJSON("/generate-paper", { topic, analysis_context: analysisContext }, 120000);

    if (data.error) {
      showResearchError(data.error);
      return;
    }

    lastResearchPaperData = data;
    renderResearchPaper(data);
  } catch (err) {
    showResearchError("Network error while generating research paper.");
  } finally {
    clearInterval(progressTimer);
    setResearchButtonLoading(false);
  }
}

function buildResearchContext() {
  let drugSimulation = null;
  try {
    const raw = localStorage.getItem("lastDrugSimulationData");
    if (raw) drugSimulation = JSON.parse(raw);
  } catch (_) {
    drugSimulation = null;
  }

  return {
    dna_analysis: lastAnalysisData || null,
    protein_analysis: lastProteinData || null,
    drug_simulation: drugSimulation || null,
  };
}

function setResearchButtonLoading(loading) {
  const btn = $("research-btn");
  if (!btn) return;

  btn.disabled = loading;
  if (loading) {
    hide("research-btn-label");
    show("research-btn-spinner");
  } else {
    show("research-btn-label");
    hide("research-btn-spinner");
  }
}

function showResearchError(msg) {
  const el = $("research-error");
  if (!el) return;
  el.textContent = "⚠ " + msg;
  show("research-error");
}

function renderResearchPaper(data) {
  const titleEl = $("research-output-title");
  const contentEl = $("research-paper-content");
  if (!titleEl || !contentEl) return;

  titleEl.textContent = `Literature Review: ${data.topic}`;
  contentEl.innerHTML = buildPaperHtml(data.paper, data.references || []);
  show("research-output");

  setTimeout(() => {
    $("research-output").scrollIntoView({ behavior: "smooth", block: "start" });
  }, 80);
}

function buildPaperHtml(rawPaper, references) {
  const sectionNames = [
    "Title",
    "Abstract",
    "Introduction",
    "Methodology",
    "Results",
    "Recent Findings",
    "Discussion",
    "Future Scope",
    "Conclusion",
    "References",
  ];

  const lines = String(rawPaper || "").split(/\r?\n/);
  const sections = {};
  let currentSection = null;

  sectionNames.forEach((name) => {
    sections[name] = [];
  });

  for (const lineRaw of lines) {
    const line = lineRaw.trim();
    if (!line) continue;

    const headingMatch = line.match(/^#{0,3}\s*\**(Title|Abstract|Introduction|Methodology|Results|Recent Findings|Discussion|Future Scope|Conclusion|References)\**\s*:?(.*)$/i);
    if (headingMatch) {
      const matchedHeading = sectionNames.find((name) => name.toLowerCase() === headingMatch[1].toLowerCase());
      currentSection = matchedHeading || null;
      const inlineText = (headingMatch[2] || "").trim();
      if (currentSection && inlineText) {
        sections[currentSection].push(inlineText);
      }
      continue;
    }

    if (currentSection) {
      sections[currentSection].push(line);
    }
  }

  let html = "";
  for (const name of sectionNames) {
    const bodyLines = sections[name];

    if (name === "References") {
      html += `<div class="research-paper-section"><h4>${escHtml(name)}</h4>`;
      if (references.length > 0) {
        html += "<ol class=\"research-ref-list\">";
        references.forEach((ref) => {
          const title = escHtml(String(ref.title || ""));
          const year = escHtml(String(ref.year || "N/A"));
          if (title) {
            html += `<li>${title} <span class=\"ref-year\">(${year})</span></li>`;
          }
        });
        html += "</ol>";
      } else {
        html += "<p>No references found.</p>";
      }
      html += "</div>";
      continue;
    }

    if (bodyLines.length > 0) {
      const text = escHtml(bodyLines.join("\n\n")).replace(/\n\n/g, "</p><p>").replace(/\n/g, "<br/>");
      html += `<div class="research-paper-section"><h4>${escHtml(name)}</h4><p>${text}</p></div>`;
    }
  }

  if (!html) {
    return `<div class="research-paper-section"><h4>Generated Paper</h4><p>${escHtml(String(rawPaper || "No content returned."))}</p></div>`;
  }

  return html;
}

async function downloadResearchPDF() {
  if (!lastResearchPaperData) {
    showResearchError("Generate a research paper first.");
    return;
  }

  const btn = $("download-research-pdf-btn");
  const label = $("research-pdf-btn-label");
  const spinner = $("research-pdf-btn-spinner");

  btn.disabled = true;
  label.classList.add("hidden");
  spinner.classList.remove("hidden");

  try {
    const res = await fetch("/api/download-research-pdf", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(lastResearchPaperData),
    });

    if (!res.ok) {
      const err = await res.json();
      showResearchError(err.error || "Failed to generate research PDF.");
      return;
    }

    const blob = await res.blob();
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `Research_Paper_${(lastResearchPaperData.topic || "topic").replace(/[^A-Za-z0-9_-]+/g, "_")}.pdf`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  } catch (err) {
    showResearchError("Network error while downloading research PDF.");
  } finally {
    btn.disabled = false;
    label.classList.remove("hidden");
    spinner.classList.add("hidden");
  }
}

// ─── Allow pressing Enter to trigger search ──────────────────
document.getElementById("gene-input").addEventListener("keydown", (e) => {
  if (e.key === "Enter") searchGene();
});

const researchInput = document.getElementById("research-input");
if (researchInput) {
  researchInput.addEventListener("keydown", (e) => {
    if (e.key === "Enter") generateResearchPaper();
  });
}
