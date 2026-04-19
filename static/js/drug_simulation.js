let comparisonChart = null;
let lastReport = "";
let dockingViewer = null;
let isSimulationRunning = false;

const pdbPool = ["1CRN", "1M17", "2ONC", "6LU7", "4HHB"];

function byId(id) {
  return document.getElementById(id);
}

function show(id) {
  byId(id).classList.remove("hidden");
}

function hide(id) {
  byId(id).classList.add("hidden");
}

function setTheme(theme) {
  document.body.className = theme === "dark" ? "theme-dark" : "theme-light";
  localStorage.setItem("dock_theme", theme);
}

function applySavedTheme() {
  const saved = localStorage.getItem("dock_theme") || "light";
  setTheme(saved);
}

function renderViewer(data) {
  const target = byId("viewer");

  if (!window.$3Dmol) {
    target.innerHTML = "<p style='color:#cbd5e1;padding:12px'>3D placeholder unavailable.</p>";
    return;
  }

  const seed = Number(data.length || 0) + Number(data.best_pose || 1);
  const pdbId = pdbPool[seed % pdbPool.length];
  if (!dockingViewer) {
    target.innerHTML = "";
    dockingViewer = $3Dmol.createViewer(target, { backgroundColor: "#020617" });
  }

  dockingViewer.clear();
  $3Dmol.download("pdb:" + pdbId, dockingViewer, { doAssembly: true }, function () {
    dockingViewer.setStyle({}, { cartoon: { color: "spectrum" } });
    dockingViewer.zoomTo();
    dockingViewer.render();
  });
}

function renderChart(data) {
  const canvas = byId("comparison-chart");
  if (comparisonChart) {
    comparisonChart.destroy();
  }

  const values = data.values || [];
  const labels = data.labels || values.map((_, idx) => idx + 1);
  const bestPose = Number(data.best_pose || 1);
  const pointColors = values.map((_, idx) => (idx + 1 === bestPose ? "#ef4444" : "#0ea5e9"));
  const pointSizes = values.map((_, idx) => (idx + 1 === bestPose ? 7 : 4));

  comparisonChart = new Chart(canvas, {
    type: "line",
    data: {
      labels,
      datasets: [
        {
          label: "Binding Energy (kcal/mol)",
          data: values,
          borderColor: "#0ea5e9",
          backgroundColor: "rgba(14, 165, 233, 0.15)",
          fill: true,
          tension: 0.28,
          pointBackgroundColor: pointColors,
          pointRadius: pointSizes,
          pointHoverRadius: 8,
        },
      ],
    },
    options: {
      responsive: true,
      plugins: {
        legend: { display: false },
        tooltip: {
          callbacks: {
            label: (ctx) => {
              const pose = ctx.dataIndex + 1;
              const suffix = pose === bestPose ? " (Best Pose)" : "";
              return ` ${ctx.raw} kcal/mol${suffix}`;
            },
          },
        },
      },
      scales: {
        x: {
          title: { display: true, text: data.x_label || "Pose Number" },
        },
        y: {
          title: { display: true, text: data.y_label || "Binding Energy (kcal/mol)" },
        },
      },
    },
  });
}

function setResult(data) {
  byId("binding-energy").textContent = data.best_energy;
  byId("avg-energy").textContent = data.average_energy;
  byId("best-pose").textContent = data.best_pose;
  byId("interaction-summary").textContent = data.interaction;
  byId("disease-name").textContent = data.target_name || data.disease || "Uploaded target";
  byId("protein-name").textContent = data.protein_name;
  byId("ligand-name").textContent = data.ligand_name;
  byId("seq-length").textContent = data.length ? data.length + " nt" : "N/A";
  byId("seq-fingerprint").textContent = data.sequence_fingerprint || "N/A";
  byId("dna-preview").textContent = data.dna_sequence_preview || "N/A";
  byId("protein-preview").textContent = data.protein_preview || "N/A";
  byId("best-frame").textContent = data.best_frame !== undefined ? data.best_frame : "N/A";
  byId("stop-position").textContent = data.stop_position !== null && data.stop_position !== undefined ? data.stop_position : "N/A";
  byId("pose-count").textContent = data.pose_count || (data.binding_energies ? data.binding_energies.length : "N/A");
  byId("ligand-type").textContent = data.ligand_type || "unknown";
  byId("best-pose-energy").textContent = data.best_energy;
  byId("avg-pose-energy").textContent = data.average_energy;
  byId("h-bonds").textContent = (data.bond_profile && data.bond_profile.hydrogen_bonds) ?? "N/A";
  byId("hydrophobic-contacts").textContent = (data.bond_profile && data.bond_profile.hydrophobic_contacts) ?? "N/A";
  byId("ionic-contacts").textContent = (data.bond_profile && data.bond_profile.ionic_interactions) ?? "N/A";
  byId("pi-stacking").textContent = (data.bond_profile && data.bond_profile.pi_stacking) ?? "N/A";
  byId("binding-pocket").textContent = (data.bond_profile && data.bond_profile.predicted_binding_pocket) || "N/A";
  byId("ncbi-topic").textContent = data.extracted_topic || (data.ncbi_context && data.ncbi_context.topic) || "N/A";
  byId("ncbi-count").textContent = (data.ncbi_context && data.ncbi_context.count) || 0;
  byId("ai-summary").textContent = data.ai_summary || "AI summary unavailable.";
  byId("ai-status").textContent = data.ai_status || "unknown";
  byId("ai-error").textContent = data.ai_error || "None";

  const papersWrap = byId("ncbi-papers");
  const papers = (data.ncbi_context && data.ncbi_context.papers) || [];
  if (!papers.length) {
    papersWrap.innerHTML = "<p class='paper-item'>No PubMed references found for extracted topic.</p>";
  } else {
    papersWrap.innerHTML = papers
      .map((p) => {
        const title = p.title || "Untitled";
        const year = p.year || "N/A";
        const journal = p.journal || "Unknown journal";
        const pmid = p.pmid || "N/A";
        return `<p class='paper-item'><strong>${title}</strong><br>${journal} (${year}) · PMID: ${pmid}</p>`;
      })
      .join("");
  }
  const modeBase = data.enrichment_enabled ? "Scientific Simulation + Enrichment" : "Scientific Simulation (Fast Mode)";
  byId("mode-badge").textContent = data.used_demo ? `${modeBase} + Demo Project` : modeBase;

  // Store simulation output for cross-page research synthesis on the main app.
  localStorage.setItem("lastDrugSimulationData", JSON.stringify(data));

  lastReport = data.report || "No report available.";
  renderViewer(data);
  renderChart(data.graph || {});
}

function updateInputModeNote() {
  const dnaSelected = byId("dna-file").files.length > 0;
  const ligandSelected = byId("ligand-file").files.length > 0;
  const demoProject = byId("demo-project").value;
  const enrichmentEnabled = byId("enable-enrichment").checked;
  const dnaInput = byId("dna-file");
  const note = byId("input-mode-note");

  const enrichTag = enrichmentEnabled ? "AI+NCBI ON" : "AI+NCBI OFF";

  if (demoProject) {
    dnaInput.required = false;
    note.textContent = `Current Mode: Demo Project selected (DNA upload optional, ${enrichTag})`;
    return;
  }

  dnaInput.required = true;

  if (!dnaSelected && !ligandSelected) {
    note.textContent = `Current Mode: Waiting for DNA file (${enrichTag})`;
    return;
  }

  if (dnaSelected && ligandSelected) {
    note.textContent = `Current Mode: DNA + Ligand uploaded (multi-pose simulation ready, ${enrichTag})`;
    return;
  }

  if (dnaSelected) {
    note.textContent = `Current Mode: DNA uploaded (ligand optional, ${enrichTag})`;
    return;
  }

  note.textContent = `Current Mode: Ligand selected (DNA file still required, ${enrichTag})`;
}

async function runSimulation(event) {
  event.preventDefault();
  if (isSimulationRunning) return;

  isSimulationRunning = true;

  hide("input-state");
  hide("result-state");
  show("processing-state");

  const form = byId("docking-form");
  const formData = new FormData(form);
  const submitBtn = form.querySelector('button[type="submit"]');
  if (submitBtn) submitBtn.disabled = true;

  try {
    const response = await fetch("/api/drug-simulation", {
      method: "POST",
      body: formData,
    });

    const data = await response.json();
    if (!response.ok) {
      throw new Error(data.error || "Simulation failed");
    }

    setTimeout(() => {
      hide("processing-state");
      show("result-state");
      setResult(data);
      byId("result-state").scrollIntoView({ behavior: "smooth", block: "start" });
    }, 250);
  } catch (error) {
    hide("processing-state");
    show("input-state");
    alert(error.message || "Failed to run simulation");
  } finally {
    isSimulationRunning = false;
    if (submitBtn) submitBtn.disabled = false;
  }
}

function downloadReport() {
  const content = lastReport || "No report available.";
  const blob = new Blob([content], { type: "text/plain" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = "Drug_Discovery_Simulation_Report.txt";
  a.click();
  URL.revokeObjectURL(url);
}

function init() {
  applySavedTheme();

  byId("start-btn").addEventListener("click", () => {
    hide("home-state");
    show("input-state");
    byId("input-state").scrollIntoView({ behavior: "smooth", block: "start" });
  });

  byId("theme-toggle").addEventListener("click", () => {
    const isDark = document.body.classList.contains("theme-dark");
    setTheme(isDark ? "light" : "dark");
  });

  byId("docking-form").addEventListener("submit", runSimulation);
  byId("dna-file").addEventListener("change", updateInputModeNote);
  byId("ligand-file").addEventListener("change", updateInputModeNote);
  byId("demo-project").addEventListener("change", updateInputModeNote);
  byId("enable-enrichment").addEventListener("change", updateInputModeNote);
  byId("download-report").addEventListener("click", downloadReport);
  byId("run-again").addEventListener("click", () => {
    hide("result-state");
    show("input-state");
    byId("docking-form").reset();
    byId("enable-enrichment").checked = true;
    updateInputModeNote();
  });

  updateInputModeNote();
}

document.addEventListener("DOMContentLoaded", init);
