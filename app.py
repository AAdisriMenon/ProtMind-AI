import streamlit as st
from backend import apply_custom_css

st.set_page_config(page_title="ProtMind AI Dashboard", page_icon="🧠", layout="centered")
apply_custom_css()

# Initialize Session State Memory to pass data between pages
if 'payload' not in st.session_state:
    st.session_state['payload'] = None

st.markdown("<h1>🧠 ProtMind AI</h1>", unsafe_allow_html=True)
st.markdown("<h3>Explainable Protein Mutation Intelligence & Therapeutic Discovery</h3>", unsafe_allow_html=True)
st.divider()

st.info("👋 Welcome to the ProtMind AI Pipeline. Please use the sidebar to navigate through the modules, starting with Step 1: User Input.")

if st.session_state['payload']:
    st.success(f"✅ Active Session Found: Analyzing {st.session_state['payload']['uniprot_id']} ({st.session_state['payload']['mutation']})")