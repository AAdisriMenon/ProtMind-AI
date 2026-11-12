import streamlit as st
import streamlit.components.v1 as components
import py3Dmol
from backend import (
    apply_custom_css,
    render_progress,
    fetch_alphafold_url,
    fetch_pdb_content,
    fetch_uniprot_data,
)

st.set_page_config(
    page_title="Step 02 | ProtMind AI", page_icon="🔍", layout="wide"
)
apply_custom_css()
render_progress(2)

st.markdown("<h2>🔍 Step 02: Data Retrieval (UniProtKB & AlphaFold DB)</h2>", unsafe_allow_html=True)

# Guard clause for missing payload
if "payload" not in st.session_state or not st.session_state["payload"]:
    st.warning(
        "⚠️ No active protein data found. Please complete Step 01: User Input & Validation first."
    )
    st.page_link("pages/1_User_Input.py", label="⬅️ Back to Step 01", icon="⚡")
    st.stop()

payload = st.session_state["payload"]
uniprot_id = payload.get("uniprot_id", "").strip()
mutation = payload.get("mutation", "").strip()

try:
    mut_pos = int(mutation[1:-1])
except ValueError:
    mut_pos = None

# Fetch & Cache live data in session_state to avoid redundant network calls on tab clicks
if "retrieved_data" not in st.session_state or st.session_state.get("current_query_id") != uniprot_id:
    with st.spinner(f"🔬 Querying UniProtKB and AlphaFold DB for {uniprot_id}..."):
        metadata = fetch_uniprot_data(uniprot_id)
        pdb_url = fetch_alphafold_url(uniprot_id)
        pdb_content = fetch_pdb_content(pdb_url) if pdb_url else ""

        st.session_state["retrieved_data"] = {
            "metadata": metadata,
            "pdb_url": pdb_url,
            "pdb_content": pdb_content,
        }
        st.session_state["current_query_id"] = uniprot_id

        # Auto-fill the sequence from UniProt if the user didn't paste/upload one in Step 01
        if not st.session_state["payload"].get("sequence") and metadata.get("Sequence"):
            st.session_state["payload"]["sequence"] = metadata["Sequence"]

retrieved = st.session_state["retrieved_data"]
uniprot_metadata = retrieved["metadata"]
af_pdb_url = retrieved["pdb_url"]
pdb_content = retrieved.get("pdb_content", "")

if payload.get("sequence"):
    st.caption(
        f"📄 Active sequence source: {payload.get('source', 'Unknown')} — "
        f"{len(payload['sequence'])} aa | Target Mutation: **{mutation}**"
    )

# ================= TABS =================
tab1, tab2 = st.tabs(["🔵 UniProtKB Annotations", "🧬 AlphaFold DB 3D Structure"])

with tab1:
    st.markdown("<h3>UniProtKB Protein Annotations & Summary</h3>", unsafe_allow_html=True)
    if "Error" in uniprot_metadata:
        st.error(uniprot_metadata["Error"])
    else:
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Protein Name", uniprot_metadata.get("Protein Name", "N/A"))
            st.metric("Gene Symbol", uniprot_metadata.get("Gene", "N/A"))
        with col2:
            st.metric("Organism", uniprot_metadata.get("Organism", "N/A"))
            st.metric("Sequence Length", f"{uniprot_metadata.get('Sequence Length', 0)} aa")

        st.markdown("#### 📖 Functional Summary")
        st.info(uniprot_metadata.get("Function", "No functional summary available."))

        col_loc, col_kw = st.columns(2)
        with col_loc:
            st.markdown("#### 📍 Subcellular Location")
            st.write(uniprot_metadata.get("Subcellular Location", "Not specified"))
        with col_kw:
            st.markdown("#### 🏷️ UniProt Keywords")
            kws = uniprot_metadata.get("Keywords", [])
            st.write(", ".join(kws) if kws else "None listed")

        with st.expander("📄 View UniProtKB Metadata JSON"):
            display_metadata = {
                k: v for k, v in uniprot_metadata.items() if k not in {"Sequence", "Features", "Pfam Cross-Refs"}
            }
            st.json(display_metadata)

with tab2:
    st.markdown("<h3>AlphaFold DB 3D Structure Visualization</h3>", unsafe_allow_html=True)
    if not pdb_content:
        if af_pdb_url:
            st.warning("⚠️ AlphaFold PDB URL located, but structure content could not be downloaded.")
            st.link_button("📥 Open .PDB URL", af_pdb_url)
        else:
            st.warning("⚠️ No AlphaFold 3D model found for this UniProt accession.")
    else:
        col_view, col_ctrl = st.columns([2, 1])

        with col_ctrl:
            st.markdown("#### 🎨 Visualization Controls")
            view_style = st.selectbox(
                "Protein Style",
                ["Cartoon (Ribbon)", "Stick", "Sphere (Space-filling)", "Line"],
            )
            color_scheme = st.selectbox(
                "Color Scheme",
                ["Default (Cyan)", "By Secondary Structure", "Spectrum (Rainbow)"],
            )
            st.info(
                f"Target residue **{mutation}** is highlighted in **RED**."
            )
            st.download_button(
                label="📥 Download AlphaFold .PDB",
                data=pdb_content,
                file_name=f"{uniprot_id}_AF.pdb",
                mime="text/plain",
            )

        with col_view:
            view = py3Dmol.view(width=700, height=500)
            view.addModel(pdb_content, "pdb")
            view.setBackgroundColor("#161A23")

            if color_scheme == "By Secondary Structure":
                color_config = {"colorscheme": "ssPyMOL"}
            elif color_scheme == "Spectrum (Rainbow)":
                color_config = {"colorscheme": "amino"}
            else:
                color_config = {"color": "#00F5D4"}

            if view_style == "Cartoon (Ribbon)":
                view.setStyle({"cartoon": color_config})
            elif view_style == "Stick":
                view.setStyle({"stick": color_config})
            elif view_style == "Sphere (Space-filling)":
                view.setStyle({"sphere": color_config})
            elif view_style == "Line":
                view.setStyle({"line": color_config})

            if mut_pos:
                view.addStyle({"resi": [mut_pos]}, {"stick": {"colorscheme": "redCarbon", "radius": 0.25}})
                view.addStyle({"resi": [mut_pos]}, {"sphere": {"color": "#FF0055", "radius": 1.2}})
                view.addLabel(
                    f"Mut: {mutation}",
                    {"fontOpacity": 1, "fontSize": 14, "fontColor": "white", "backgroundColor": "#FF0055"},
                    {"resi": [mut_pos]},
                )
                view.zoomTo({"resi": [mut_pos]})
            else:
                view.zoomTo()

            html_string = f'<div style="height: 500px; width: 100%;">{view._make_html()}</div>'
            components.html(html_string, height=520)

st.markdown("---")
if not payload.get("sequence"):
    st.error(
        "❌ No sequence could be retrieved automatically. Please go back to Step 01 and "
        "paste or upload a sequence manually."
    )
    st.page_link("pages/1_User_Input.py", label="⬅️ Back to Step 01", icon="⚡")
else:
    st.page_link(
        "pages/3_Preprocessing.py",
        label="Proceed to Step 03: Preprocessing & Motif/Domain Search ➡️",
        icon="⚙️",
    )
