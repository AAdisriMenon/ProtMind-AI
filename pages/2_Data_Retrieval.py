import streamlit as st
import streamlit.components.v1 as components
import py3Dmol
from backend import (
    apply_custom_css,
    render_progress,
    fetch_alphafold_prediction,
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
except (ValueError, TypeError):
    mut_pos = None

# Fetch & Cache live data in session_state to avoid redundant network calls on tab clicks
if "retrieved_data" not in st.session_state or st.session_state.get("current_query_id") != uniprot_id:
    with st.spinner(f"🔬 Querying UniProtKB and AlphaFold DB for {uniprot_id}..."):
        metadata = fetch_uniprot_data(uniprot_id)
        af_pred = fetch_alphafold_prediction(uniprot_id, target_pos=mut_pos)
        pdb_url = af_pred.get("pdbUrl") if af_pred else fetch_alphafold_url(uniprot_id)
        pdb_content = fetch_pdb_content(pdb_url) if pdb_url else ""

        st.session_state["retrieved_data"] = {
            "metadata": metadata,
            "af_pred": af_pred,
            "pdb_url": pdb_url,
            "pdb_content": pdb_content,
        }
        st.session_state["current_query_id"] = uniprot_id

        # Auto-fill the sequence from UniProt if the user didn't paste/upload one in Step 01
        if not st.session_state["payload"].get("sequence") and metadata.get("Sequence"):
            st.session_state["payload"]["sequence"] = metadata["Sequence"]

retrieved = st.session_state["retrieved_data"]
uniprot_metadata = retrieved["metadata"]
af_pred = retrieved.get("af_pred", {})
af_pdb_url = retrieved.get("pdb_url", "")
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
    st.markdown("<h3>AlphaFold DB 3D Structure (Entire Sequence)</h3>", unsafe_allow_html=True)
    if not pdb_content:
        if af_pdb_url:
            st.warning("⚠️ AlphaFold PDB URL located, but structure content could not be downloaded.")
            st.link_button("📥 Open .PDB URL", af_pdb_url)
        else:
            st.warning("⚠️ No AlphaFold 3D model found for this UniProt accession.")
    else:
        # Check if model has multiple fragments (for large proteins)
        all_frags = af_pred.get("allFragments", []) if isinstance(af_pred, dict) else []
        if len(all_frags) > 1:
            frag_labels = [
                f"{f.get('entryId', f'Fragment {i+1}')} (Residues {f.get('uniprotStart', 1)}–{f.get('uniprotEnd', '?')})"
                for i, f in enumerate(all_frags)
            ]
            default_idx = 0
            if mut_pos:
                for i, f in enumerate(all_frags):
                    if f.get("uniprotStart", 0) <= mut_pos <= f.get("uniprotEnd", 0):
                        default_idx = i
                        break
            selected_frag_label = st.selectbox("🌐 Select AlphaFold Structure Segment:", frag_labels, index=default_idx)
            selected_frag_idx = frag_labels.index(selected_frag_label)
            active_frag = all_frags[selected_frag_idx]
            if active_frag.get("pdbUrl") != af_pdb_url:
                af_pdb_url = active_frag.get("pdbUrl")
                pdb_content = fetch_pdb_content(af_pdb_url)
                st.session_state["retrieved_data"]["pdb_url"] = af_pdb_url
                st.session_state["retrieved_data"]["pdb_content"] = pdb_content
                st.session_state["retrieved_data"]["af_pred"] = active_frag
                af_pred = active_frag

        # Sequence coverage and model metadata
        seq_len = len(payload.get("sequence", "")) or uniprot_metadata.get("Sequence Length", 0)
        cov_start = af_pred.get("uniprotStart", 1) if isinstance(af_pred, dict) else 1
        cov_end = af_pred.get("uniprotEnd", seq_len) if isinstance(af_pred, dict) else seq_len
        model_id = af_pred.get("entryId", f"AF-{uniprot_id}-F1") if isinstance(af_pred, dict) else f"AF-{uniprot_id}-F1"
        is_entire = (cov_start == 1 and (cov_end == seq_len or seq_len == 0))

        m_col1, m_col2, m_col3, m_col4 = st.columns(4)
        with m_col1:
            st.metric("AlphaFold Model", model_id)
        with m_col2:
            st.metric("Structure Scope", "Entire Sequence" if is_entire else f"Segment {cov_start}–{cov_end}")
        with m_col3:
            st.metric("Residues Covered", f"{cov_start} – {cov_end} ({cov_end - cov_start + 1} aa)")
        with m_col4:
            coverage_pct = round(((cov_end - cov_start + 1) / seq_len) * 100, 1) if seq_len else 100
            st.metric("Sequence Coverage", f"{coverage_pct}%")

        col_view, col_ctrl = st.columns([2, 1])

        with col_ctrl:
            st.markdown("#### 🎨 Visualization Controls")

            camera_options = ["Entire Sequence (Full Structure)"]
            if mut_pos:
                camera_options.append(f"Focus on Mutation Site ({mutation})")

            camera_focus = st.radio(
                "📷 Camera Focus",
                camera_options,
                index=0,
                help="Select 'Entire Sequence' to frame the full 3D folded protein chain, or 'Focus on Mutation Site' to zoom in tightly to the mutation."
            )

            view_style = st.selectbox(
                "Protein Representation",
                ["Cartoon (Ribbon)", "Stick", "Sphere (Space-filling)", "Line"],
            )
            color_scheme = st.selectbox(
                "Color Scheme",
                [
                    "AlphaFold Confidence (pLDDT)",
                    "Spectrum / Rainbow (N ➔ C)",
                    "By Secondary Structure",
                    "ProtMind Cyan",
                ],
            )

            highlight_mutation = st.checkbox("Highlight Target Residue", value=True)
            if mut_pos:
                st.caption(f"📍 Target residue **{mutation}** (pos {mut_pos}) is highlighted in **RED**.")

            st.download_button(
                label="📥 Download AlphaFold .PDB",
                data=pdb_content,
                file_name=f"{uniprot_id}_entire_sequence_AF.pdb",
                mime="text/plain",
            )

            if color_scheme == "AlphaFold Confidence (pLDDT)":
                st.markdown("""
                **AlphaFold pLDDT Legend:**
                - 🟦 **> 90**: Very high confidence
                - 🩵 **70 – 90**: Confident
                - 🟨 **50 – 70**: Low confidence
                - 🟧 **< 50**: Very low confidence (disordered)
                """)

        with col_view:
            view = py3Dmol.view(width=720, height=520)
            view.addModel(pdb_content, "pdb")
            view.setBackgroundColor("#161A23")

            if color_scheme == "AlphaFold Confidence (pLDDT)":
                color_config = {"colorscheme": {"prop": "b", "gradient": "roygb", "min": 50, "max": 90}}
            elif color_scheme == "By Secondary Structure":
                color_config = {"colorscheme": "ssPyMOL"}
            elif color_scheme == "Spectrum / Rainbow (N ➔ C)":
                color_config = {"colorscheme": "amino"}
            else:
                color_config = {"color": "#00F5D4"}

            if view_style == "Cartoon (Ribbon)":
                view.setStyle({}, {"cartoon": color_config})
            elif view_style == "Stick":
                view.setStyle({}, {"stick": color_config})
            elif view_style == "Sphere (Space-filling)":
                view.setStyle({}, {"sphere": color_config})
            elif view_style == "Line":
                view.setStyle({}, {"line": color_config})

            if mut_pos and highlight_mutation:
                view.addStyle({"resi": [mut_pos]}, {"stick": {"colorscheme": "redCarbon", "radius": 0.28}})
                view.addStyle({"resi": [mut_pos]}, {"sphere": {"color": "#FF0055", "radius": 1.2}})
                view.addLabel(
                    f"Mut: {mutation}",
                    {"fontOpacity": 1, "fontSize": 13, "fontColor": "white", "backgroundColor": "#FF0055"},
                    {"resi": [mut_pos]},
                )

            # Camera zoom: defaults to Entire Sequence (view.zoomTo()) so whole protein is visible
            if camera_focus.startswith("Focus on Mutation") and mut_pos:
                view.zoomTo({"resi": [mut_pos]})
            else:
                view.zoomTo()

            html_string = f'<div style="height: 520px; width: 100%; border-radius: 12px; overflow: hidden;">{view._make_html()}</div>'
            components.html(html_string, height=540)

        # Full sequence view expander
        if payload.get("sequence"):
            with st.expander(f"📄 View Entire Amino Acid Sequence ({len(payload['sequence'])} aa)", expanded=False):
                full_seq = payload["sequence"]
                formatted_seq = "\n".join([full_seq[i:i+60] for i in range(0, len(full_seq), 60)])
                st.code(formatted_seq, language="text")

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
