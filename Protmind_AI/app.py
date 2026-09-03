import streamlit as st
from backend import apply_custom_css, PIPELINE_STEPS

st.set_page_config(page_title="ProtMind AI Dashboard", page_icon="🧠", layout="centered")
apply_custom_css()

# Initialize Session State Memory to pass data between pages
if 'payload' not in st.session_state:
    st.session_state['payload'] = None

st.markdown("<h1>🧠 ProtMind AI</h1>", unsafe_allow_html=True)
st.markdown("<h3>Explainable Protein Mutation Intelligence & Therapeutic Discovery</h3>", unsafe_allow_html=True)
st.divider()

st.info("👋 Welcome to the ProtMind AI Pipeline. Please use the sidebar to navigate through the modules, starting with Step 1: User Input.")

with st.expander("🧭 What does this pipeline do?", expanded=False):
    for i, step in enumerate(PIPELINE_STEPS, start=1):
        st.markdown(f"**{i}. {step}**")
    st.caption(
        "Enter a UniProt ID and a point mutation (e.g. P04637 / R175H) and ProtMind AI will "
        "retrieve the protein record, compute biophysical feature shifts, estimate pathogenicity, "
        "pull known interaction partners, and render the mutated residue on the AlphaFold structure."
    )

if st.session_state['payload']:
    st.success(
        f"✅ Active Session Found: Analyzing {st.session_state['payload']['uniprot_id']} "
        f"({st.session_state['payload']['mutation']})"
    )
    st.page_link("pages/2_Data_Retrieval.py", label="Continue to Step 2 - Data Retrieval", icon="🔍")
else:
    st.page_link("pages/1_User_Input.py", label="🚀 Start Analysis: Step 1 - User Input", icon="⚡")
