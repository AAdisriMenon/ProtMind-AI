import streamlit as st
from backend import apply_custom_css

st.set_page_config(page_title="ProtMind AI", page_icon="🧬", layout="centered")
apply_custom_css()

# Initialize Session State Memory to pass data between pages
if 'payload' not in st.session_state:
    st.session_state['payload'] = None

st.markdown("<h1>🧬 ProtMind AI</h1>", unsafe_allow_html=True)
st.markdown("<h3>Protein Structural & Functional Genomics Workbench</h3>", unsafe_allow_html=True)
st.divider()

st.info("👋 Welcome to **ProtMind AI**. Follow the 3-step pipeline below to analyze your protein variant:")

with st.expander("🧭 3-Step Workflow Architecture", expanded=True):
    st.markdown("""
- **Step 01 – User Input & Validation**: Provide UniProt ID, mutation (e.g., R175H), and protein sequence (FASTA/upload).
- **Step 02 – Data Retrieval**: Fetch protein annotations and 3D structure from UniProtKB and AlphaFold DB.
- **Step 03 – Preprocessing & Motif/Domain Search**: Extract sequence features, identify motifs and conserved domains (Pfam).
    """)

if st.session_state['payload']:
    st.success(
        f"✅ Active Session Found: Analyzing {st.session_state['payload']['uniprot_id']} "
        f"({st.session_state['payload']['mutation']})"
    )
    st.page_link("pages/2_Data_Retrieval.py", label="Continue to Step 02 – Data Retrieval", icon="🔍")
else:
    st.page_link("pages/1_User_Input.py", label="🚀 Start Analysis: Step 01 – User Input & Validation", icon="⚡")
