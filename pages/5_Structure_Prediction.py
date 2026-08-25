import streamlit as st
import time
import py3Dmol
from stmol import showmol
from backend import apply_custom_css, fetch_alphafold_url, fetch_pdb_content

st.set_page_config(page_title="Step 5 | ProtMind AI", page_icon="🧩")
apply_custom_css()

st.markdown("<h2>🧩 Step 5: Protein Structure Prediction</h2>", unsafe_allow_html=True)

if not st.session_state.get('payload'):
    st.warning("⚠️ No active protein data found. Please complete Step 1: User Input first.")
    st.stop()

payload = st.session_state['payload']
uniprot_id = payload['uniprot_id']
mutation = payload['mutation']

# Extract the mutation position number for 3D highlighting
try:
    mut_pos = int(mutation[1:-1])
except ValueError:
    mut_pos = None

with st.spinner(f"🔬 Fetching AlphaFold 3D structure and biophysical data for {uniprot_id}..."):
    af_pdb_url = fetch_alphafold_url(uniprot_id)
    pdb_content = fetch_pdb_content(af_pdb_url)
    time.sleep(1)

if not pdb_content:
    st.error(f"❌ Could not retrieve 3D structure for {uniprot_id}. It may not be available in the AlphaFold database.")
else:
    col1, col2 = st.columns([2, 1])

    with col1:
        st.markdown("<h3>🔍 Interactive 3D Viewer</h3>", unsafe_allow_html=True)
        
        # Initialize py3Dmol Viewer
        view = py3Dmol.view(width=500, height=500)
        view.addModel(pdb_content, 'pdb')
        
        # Set background to match your premium dark UI
        view.setBackgroundColor('#161A23')
        
        # Style the main protein backbone
        view.setStyle({'cartoon': {'color': '#00F5D4'}})
        
        # Highlight the specific mutated residue
        if mut_pos:
            view.addStyle({'resi': str(mut_pos)}, {'stick': {'colorscheme': 'redCarbon', 'radius': 0.2}})
            view.addStyle({'resi': str(mut_pos)}, {'sphere': {'color': '#FF0055', 'radius': 1.5}})
            view.addLabel(
                f"Mut: {mutation}", 
                {'fontOpacity': 1, 'fontSize': 14, 'fontColor': 'white', 'backgroundColor': '#FF0055'}, 
                {'resi': str(mut_pos)}
            )
            # Auto-zoom directly to the mutation site
            view.zoomTo({'resi': str(mut_pos)})
        else:
            view.zoomTo()
            
        # Render the py3Dmol object in Streamlit
        showmol(view, height=500, width=500)

    with col2:
        st.markdown("<h3>📊 Structural Insights</h3>", unsafe_allow_html=True)
        st.info("The 3D conformation is generated from AlphaFold. The target residue is highlighted in **RED**.")
        
        st.metric(label="Structure Source", value="AlphaFold2")
        st.metric(label="Target Mutation", value=mutation)
        
        if mut_pos:
            st.write(f"**Visualization:** Camera is locked onto position **{mut_pos}** to evaluate steric clashes and pocket disruptions.")
        
        st.markdown("<br>", unsafe_allow_html=True)
        st.download_button(
            label="📥 Download .PDB File",
            data=pdb_content,
            file_name=f"{uniprot_id}_AF.pdb",
            mime="text/plain"
        )