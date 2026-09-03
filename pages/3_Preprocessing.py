import streamlit as st
from backend import (
    apply_custom_css,
    render_progress,
    calculate_alignment_score,
    extract_features,
)

st.set_page_config(
    page_title="Step 3 | ProtMind AI", page_icon="⚙️", layout="wide"
)
apply_custom_css()
render_progress(3)

st.markdown(
    "<h2>⚙️ Step 3: Preprocessing & Feature Extraction</h2>",
    unsafe_allow_html=True,
)

# Guard clause for payload
if "payload" not in st.session_state or not st.session_state["payload"]:
    st.warning(
        "⚠️ No active protein data found. Please complete Step 1: User Input first."
    )
    st.page_link("pages/1_User_Input.py", label="⬅️ Back to Step 1", icon="⚡")
    st.stop()

payload = st.session_state["payload"]
raw_seq = payload.get("sequence", "").strip()
mutation = payload.get("mutation", "").strip()

if not raw_seq:
    st.error(
        "❌ Sequence is missing. Please ensure sequence retrieval or entry succeeded in previous steps."
    )
    st.page_link("pages/2_Data_Retrieval.py", label="⬅️ Back to Step 2", icon="🔍")
    st.stop()

# Cache preprocessing computations in session_state
if "preprocessed_data" not in st.session_state or st.session_state.get(
    "processed_mutation"
) != mutation:
    with st.spinner("⏳ Extracting physicochemical features and computing alignment..."):
        try:
            features = extract_features(raw_seq, mutation)
        except ValueError:
            st.error(f"❌ Could not parse mutation position from '{mutation}'. Expected format like R175H.")
            st.stop()

        # Build mutant sequence only if position and wildtype match
        if features.get("Positional Match"):
            mutant_list = list(raw_seq)
            mutant_list[features["Position"] - 1] = features["Mutant AA"]
            mutant_sequence = "".join(mutant_list)
        else:
            mutant_sequence = raw_seq

        alignment_score = calculate_alignment_score(raw_seq, mutant_sequence)

        st.session_state["preprocessed_data"] = {
            "features": features,
            "mutant_sequence": mutant_sequence,
            "alignment_score": alignment_score,
        }
        st.session_state["processed_mutation"] = mutation

proc_data = st.session_state["preprocessed_data"]
features = proc_data["features"]
alignment_score = proc_data["alignment_score"]

# ================= DISPLAY COLUMNS =================
col_feat, col_align = st.columns(2)

with col_feat:
    st.markdown("<h3>Physicochemical Analysis</h3>", unsafe_allow_html=True)

    if not features.get("In Range"):
        st.error(
            f"⚠️ Mutation position {features.get('Position')} is outside the sequence "
            f"(length {len(raw_seq)}). Double check the mutation and sequence match."
        )
    elif not features.get("Positional Match"):
        st.error(
            f"⚠️ Positional mismatch: Expected `{features.get('Wildtype AA')}` at position {features.get('Position')}, but found a different residue."
        )
    else:
        st.success("✅ Mutation position verified against wildtype sequence.")

    f_col1, f_col2 = st.columns(2)
    with f_col1:
        st.metric("Wildtype Residue", f"{features['Wildtype AA']} ({features['WT Property']})")
    with f_col2:
        st.metric("Mutant Residue", f"{features['Mutant AA']} ({features['Mutant Property']})")

    with st.expander("📄 View Full Feature Dictionary"):
        st.json(features)

with col_align:
    st.markdown("<h3>Sequence & Alignment Metrics</h3>", unsafe_allow_html=True)

    st.metric(label="Sequence Identity Score", value=f"{alignment_score:.2f}%")
    st.progress(min(max(alignment_score / 100.0, 0.0), 1.0))

    if features["WT Property"] != features["Mutant Property"]:
        st.warning(f"**Biophysical Shift:** {features['Property Shift']}")
    else:
        st.info(f"**Conserved Class:** Maintained as {features['WT Property']}")

st.markdown("---")
# Fixed: this used to point to a page ("pages/4_Prediction.py") that doesn't
# exist in this project, so the button silently fell through every time.
st.page_link("pages/4_Mutation_Analysis.py", label="Proceed to Step 4: Mutation Analysis ➡️", icon="🧬")
