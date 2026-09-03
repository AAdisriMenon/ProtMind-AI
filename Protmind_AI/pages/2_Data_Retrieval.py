import streamlit as st
from backend import apply_custom_css, render_progress, fetch_alphafold_url, fetch_uniprot_data

st.set_page_config(
    page_title="Step 2 | ProtMind AI", page_icon="🔍", layout="wide"
)
apply_custom_css()
render_progress(2)

st.markdown("<h2>🔍 Step 2: Data Retrieval Results</h2>", unsafe_allow_html=True)

# Guard clause for missing payload
if "payload" not in st.session_state or not st.session_state["payload"]:
    st.warning(
        "⚠️ No active protein data found. Please complete Step 1: User Input first."
    )
    st.page_link("pages/1_User_Input.py", label="⬅️ Back to Step 1", icon="⚡")
    st.stop()

payload = st.session_state["payload"]
uniprot_id = payload.get("uniprot_id", "").strip()

# Fetch & Cache live data in session_state to avoid redundant network calls on tab clicks
if "retrieved_data" not in st.session_state or st.session_state.get(
    "current_query_id"
) != uniprot_id:
    with st.spinner(f"🔬 Querying live databases for {uniprot_id}..."):
        metadata = fetch_uniprot_data(uniprot_id)
        pdb_url = fetch_alphafold_url(uniprot_id)

        st.session_state["retrieved_data"] = {
            "metadata": metadata,
            "pdb_url": pdb_url,
        }
        st.session_state["current_query_id"] = uniprot_id

        # Auto-fill the sequence from UniProt if the user didn't paste/upload one in Step 1
        if not st.session_state["payload"].get("sequence") and metadata.get("Sequence"):
            st.session_state["payload"]["sequence"] = metadata["Sequence"]

retrieved = st.session_state["retrieved_data"]
uniprot_metadata = retrieved["metadata"]
af_pdb_url = retrieved["pdb_url"]

if payload.get("sequence"):
    st.caption(f"📄 Active sequence source: {payload.get('source', 'Unknown')} — {len(payload['sequence'])} aa")

# ================= TABS =================
tab1, tab2, tab3, tab4 = st.tabs(
    ["🔵 UniProt", "🏥 NIH / ClinVar", "👥 gnomAD", "🧬 Structure"]
)

with tab1:
    st.markdown("<h3>Protein Knowledgebase</h3>", unsafe_allow_html=True)
    if "Error" in uniprot_metadata:
        st.error(uniprot_metadata["Error"])
    else:
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Protein Name", uniprot_metadata.get("Protein Name", "N/A"))
            st.metric("Gene Symbol", uniprot_metadata.get("Gene", "N/A"))
        with col2:
            st.metric("Organism", uniprot_metadata.get("Organism", "N/A"))
            st.metric(
                "Sequence Length",
                f"{uniprot_metadata.get('Sequence Length', 0)} aa",
            )

        with st.expander("📄 View Raw Metadata JSON"):
            # Sequence text can be long; keep the quick-view JSON compact
            display_metadata = {k: v for k, v in uniprot_metadata.items() if k != "Sequence"}
            st.json(display_metadata)

with tab2:
    st.markdown("<h3>Clinical Significance</h3>", unsafe_allow_html=True)
    col1, col2 = st.columns(2)
    with col1:
        st.metric(
            label="Target Mutation", value=payload.get("mutation", "None")
        )
    with col2:
        st.metric(label="ClinVar Status", value="Under Review")
    st.caption("ℹ️ Live ClinVar lookup is not wired up yet — this is a placeholder for the next integration.")

with tab3:
    st.markdown("<h3>Population Genomics</h3>", unsafe_allow_html=True)
    allele_frequencies = {
        "European": 0.005,
        "African": 0.001,
        "Asian": 0.003,
        "Global": 0.004,
    }
    st.bar_chart(allele_frequencies)
    st.caption("ℹ️ Sample/illustrative allele frequencies — swap in a live gnomAD API call for production.")

with tab4:
    st.markdown("<h3>Predicted 3D Structure</h3>", unsafe_allow_html=True)
    if af_pdb_url:
        st.success("✅ AlphaFold structure located!")
        st.link_button("📥 Download .PDB File", af_pdb_url)
    else:
        st.warning("⚠️ No AlphaFold 3D model found for this accession.")

st.markdown("---")
if not payload.get("sequence"):
    st.error(
        "❌ No sequence could be retrieved automatically. Please go back to Step 1 and "
        "paste or upload a sequence manually."
    )
    st.page_link("pages/1_User_Input.py", label="⬅️ Back to Step 1", icon="⚡")
else:
    st.page_link("pages/3_Preprocessing.py", label="Proceed to Step 3: Preprocessing ➡️", icon="⚙️")
