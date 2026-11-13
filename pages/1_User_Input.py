import streamlit as st
from backend import apply_custom_css, render_progress, validate_sequence, validate_uniprot_id, validate_mutation, parse_fasta

st.set_page_config(page_title="Step 01 | ProtMind AI", page_icon="⚡")
apply_custom_css()
render_progress(1)

st.markdown("<h2>⚡ Step 01: User Input & Validation</h2>", unsafe_allow_html=True)

with st.form("user_input_form", clear_on_submit=False):
    col1, col2 = st.columns(2)
    with col1:
        uniprot_id = st.text_input("🧬 UniProt ID", placeholder="P04637", help="Example: P04637")
    with col2:
        mutation = st.text_input("🧪 Mutation", placeholder="R175H", help="Format: [WT][Position][Mut] (e.g., R175H)")

    st.markdown("---")
    st.markdown("### 📂 Sequence Input Options (optional)")
    st.caption(
        "You can leave this blank — ProtMind AI will automatically pull the canonical "
        "sequence from UniProt in Step 2. Provide it manually only if you want to analyze "
        "a custom/non-canonical sequence."
    )

    # Option A: File Upload
    uploaded_file = st.file_uploader(
        "Upload FASTA File (.fasta, .fa, .txt)",
        type=["fasta", "fa", "txt"],
        help="Upload a FASTA file containing your target protein sequence."
    )

    st.markdown("<p style='text-align:center;'>— OR —</p>", unsafe_allow_html=True)

    # Option B: Manual Sequence Paste
    sequence_input = st.text_area(
        "Paste Protein Sequence or FASTA text",
        placeholder=">sp|P04637|...\nMEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPL...",
        height=150
    )

    submitted = st.form_submit_button("🚀 Validate & Save Data")

if submitted:
    errors = []
    final_sequence = ""
    source_type = ""

    # Priority 1: File Upload
    if uploaded_file is not None:
        try:
            raw_text = uploaded_file.getvalue().decode("utf-8")
            final_sequence = parse_fasta(raw_text)
            source_type = f"Uploaded file ({uploaded_file.name})"
        except Exception:
            errors.append("❌ Could not read the uploaded FASTA file. Ensure it is UTF-8 encoded text.")
    # Priority 2: Text Area
    elif sequence_input.strip():
        final_sequence = parse_fasta(sequence_input)
        source_type = "Manual input"
    # Priority 3: No sequence provided -- will be auto-fetched from UniProt in Step 2
    else:
        final_sequence = ""
        source_type = "Auto-fetch from UniProt (Step 2)"

    # Validation Checks
    if final_sequence and not validate_sequence(final_sequence):
        errors.append("❌ Invalid amino acid characters detected. Only standard 20 amino acids allowed.")

    if not uniprot_id or not validate_uniprot_id(uniprot_id):
        errors.append("❌ Invalid UniProt ID format. Example: P04637")

    if not mutation or not validate_mutation(mutation):
        errors.append("❌ Invalid Mutation format. Example: R175H")

    # Output Handling
    if errors:
        for error in errors:
            st.error(error)
    else:
        st.session_state['payload'] = {
            "sequence": final_sequence.strip().upper(),
            "uniprot_id": uniprot_id.strip().upper(),
            "mutation": mutation.strip().upper(),
            "source": source_type
        }
        st.success(f"✅ Input validated! Sequence source: {source_type}")
        if final_sequence:
            st.info(f"**Loaded Sequence Length:** {len(final_sequence)} amino acids")
        else:
            st.info("**Sequence will be fetched automatically from UniProt in the next step.**")
        st.page_link("pages/2_Data_Retrieval.py", label="👉 Proceed to Step 2 (Data Retrieval)", icon="🔍")
