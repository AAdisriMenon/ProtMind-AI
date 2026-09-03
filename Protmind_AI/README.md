# 🧠 ProtMind AI

AI-powered platform for protein mutation analysis and prediction using bioinformatics and machine learning.

Enter a UniProt ID and a point mutation (e.g. `P04637` / `R175H`) and ProtMind AI walks through a
5-step pipeline: input validation → live UniProt/AlphaFold retrieval → biophysical feature
extraction → pathogenicity + interaction-network analysis → interactive 3D structure viewer.

## Run locally

```bash
git clone <your-repo-url>
cd Protmind_AI
python -m venv .venv
source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
streamlit run app.py
```

Open the URL Streamlit prints (usually `http://localhost:8501`).

## Deploy (Streamlit Community Cloud)

1. Push this folder to a GitHub repo (make sure `.venv/` is **not** committed — see `.gitignore`).
2. Go to [share.streamlit.io](https://share.streamlit.io) → **New app**.
3. Pick your repo/branch, set **Main file path** to `app.py`, and deploy.
4. Streamlit Cloud reads `requirements.txt` and `.streamlit/config.toml` automatically.

## What changed in this version

- Fixed a UniProt-ID regex bug where the alternation wasn't properly anchored, so malformed IDs
  could slip past validation.
- Fixed an `IndexError` crash when a mutation's position was outside the sequence length.
- Fixed the "Sequence" auto-fill: Step 2 fetches the real UniProt sequence now (previously the key
  it looked for was never populated, so it was dead code).
- Made manual sequence entry in Step 1 **optional** — it auto-fetches from UniProt in Step 2 if omitted.
- Fixed a broken "Proceed" button in Step 3 that pointed at a page that didn't exist.
- Replaced the hardcoded pathogenicity scores in Step 4 with a real, explainable heuristic based on
  the BLOSUM62 substitution matrix + physicochemical property shifts (clearly labeled as a heuristic,
  not a clinical-grade predictor).
- Added `@st.cache_data` to all network calls so repeated navigation doesn't re-hit external APIs.
- Added a sidebar pipeline-progress tracker and working forward/back navigation on every step.
- Added a missing `pandas` dependency to `requirements.txt` (Step 4 imports it but it was never
  listed — this would have broken a fresh deploy) and removed the abandoned `ipython-genutils`
  package.
- Added `.streamlit/config.toml` so the dark theme applies immediately instead of flashing light mode.

## Notes for judges / reviewers

- Pathogenicity scores (Step 4) and gnomAD/ClinVar panels (Step 2) are clearly labeled where they
  are illustrative/heuristic rather than live clinical database calls — this keeps the demo honest
  about what's real-time (UniProt, AlphaFold, STRING) vs. simulated.
