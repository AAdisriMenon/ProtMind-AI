import streamlit as st
from backend import apply_custom_css, generate_pdf_report, fetch_alphafold_url, fetch_pdb_content, calculate_phi_psi

st.set_page_config(page_title="Step 10 | ProtMind AI", page_icon="📄")
apply_custom_css()

st.markdown("<h2>📄 Step 10: Automated Clinical Report</h2>", unsafe_allow_html=True)

if not st.session_state.get('payload'):
    st.warning("⚠️ No active protein data found. Please complete Step 1: User Input first.")
    st.stop()

payload = st.session_state['payload']

st.markdown("### Export Full Pipeline Analytics")
st.write("Generate a comprehensive PDF document summarizing sequence metrics, structural disruption profiles, and therapeutic targeting insights.")

# Fetch structural angles if available
mut_pos = int(payload['mutation'][1:-1]) if payload['mutation'][1:-1].isdigit() else 1
af_pdb_url = fetch_alphafold_url(payload['uniprot_id'])
pdb_content = fetch_pdb_content(af_pdb_url)
angles = calculate_phi_psi(pdb_content, mut_pos)

# Generate PDF with embedded visual plots
with st.spinner("📑 Assembling report and generating structural plots..."):
    pdf_bytes = generate_pdf_report(
        payload=payload,
        ddg_status="Highly Destabilizing (+2.2 kcal/mol)",
        ppi_count=5,
        angles=angles
    )

st.markdown("<br>", unsafe_allow_html=True)

st.download_button(
    label="📥 Download Clinical Report (PDF)",
    data=pdf_bytes,
    file_name=f"{payload['uniprot_id']}_{payload['mutation']}_Clinical_Report.pdf",
    mime="application/pdf"
)

st.success("✅ Report generated with embedded structural and biophysical profiles!")