# dna-analyzer2

A Flask-based DNA and protein bioinformatics app with sequence analysis, docking-style multi-pose simulation, NCBI context enrichment, and AI-generated interpretation.

## Features
- DNA analysis: GC%, length, reverse complement, ORF scan (3 frames), protein preview
- Simulated docking: 10-15 poses, best/average energy, pose chart
- Interaction profile: estimated hydrogen bonds, hydrophobic/ionic/pi-stacking contacts
- NCBI enrichment: PubMed context by extracted target topic
- AI insight: one detailed paragraph summary (with robust fallback)
- Demo projects: run without file upload (`BRCA1`, `TP53`, `EGFR`)

## Run Locally

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
python app.py
```

Open:
- Main app: http://127.0.0.1:5000
- Simulation page: http://127.0.0.1:5000/drug-discovery

## Environment Variables
- `ENTREZ_EMAIL` required for NCBI Entrez usage
- `GEMINI_API_KEY` required for AI summary generation
- `GEMINI_MODEL` optional (default fallback chain is used if omitted)

## Deploy

Use `render.yaml` with Render Blueprint deploy, or configure:
- Build: `pip install -r requirements.txt`
- Start: `gunicorn app:app`
