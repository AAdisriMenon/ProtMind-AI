# 🧬 ProtMind AI — Protein Structural & Functional Genomics Workbench

AI-powered platform for protein mutation analysis and functional genomics using bioinformatics pipelines, live database queries, and structural visualization.

Enter a UniProt ID and a point mutation (e.g. `P04637` / `R175H`) and ProtMind AI walks through an intuitive **3-step workflow**: 
1. **User Input & Validation** → 
2. **Data Retrieval (UniProtKB & AlphaFold DB)** → 
3. **Preprocessing & Motif/Domain Analysis (Pfam & Biophysical Profiling)**.

---

## 🚀 Run Locally

```bash
git clone [https://github.com/AAdisriMenon/ProtMind-AI.git](https://github.com/AAdisriMenon/ProtMind-AI.git)
cd ProtMind-AI
python -m venv venv
venv\Scripts\activate      # On Windows PowerShell
pip install -r requirements.txt
streamlit run app.py
