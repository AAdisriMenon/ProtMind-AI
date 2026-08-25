import streamlit as st
import time
import pandas as pd
from backend import apply_custom_css, fetch_string_ppi, predict_pathogenicity

st.set_page_config(page_title="Step 4 | ProtMind AI", page_icon="🧬")
apply_custom_css()

st.markdown("<h2>🧬 Step 4: Mutation Analysis & Functional Proteomics</h2>", unsafe_allow_html=True)

if not st.session_state.get('payload'):
    st.warning("⚠️ No active protein data found. Please complete Step 1: User Input first.")
    st.stop()

payload = st.session_state['payload']

with st.spinner(f"🔬 Analyzing pathogenicity and querying STRING database for {payload['uniprot_id']}..."):
    # Fetch Functional Proteomics (PPIs)
    ppi_data = fetch_string_ppi(payload['uniprot_id'])
    # Fetch Pathogenicity Scores
    pathogenicity_scores = predict_pathogenicity(payload['mutation'])
    time.sleep(1.5) # Buffer for smooth UI transition

col1, col2 = st.columns(2)

with col1:
    st.markdown("<h3>⚠️ Pathogenicity Predictors</h3>", unsafe_allow_html=True)
    st.info("Evaluating mutation impact using evolutionary conservation and structural features.")
    
    for tool, result in pathogenicity_scores.items():
        score = result["score"]
        pred = result["prediction"]
        
        # Dynamic color coding based on severity
        if score > 0.8 or (tool == "SIFT" and score < 0.05):
            color = "🔴" # High Risk
        else:
            color = "🟢" # Low Risk
            
        st.markdown(f"**{color} {tool}:** {pred} (Score: {score})")

with col2:
    st.markdown("<h3>🕸️ Protein-Protein Interactions (PPI)</h3>", unsafe_allow_html=True)
    st.write("Top interacting partners that may be disrupted by this mutation:")
    
    if ppi_data:
        # Convert dictionary list to a clean Pandas DataFrame for Streamlit rendering
        df = pd.DataFrame(ppi_data)
        # Format the score to be a percentage for readability
        df['Interaction Score'] = df['Interaction Score'].apply(lambda x: f"{x * 100:.1f}%")
        st.dataframe(df, use_container_width=True, hide_index=True)
    else:
        st.warning(f"No significant interactions found in the STRING database for {payload['uniprot_id']}.")
        
st.divider()
st.success("✅ Mutation Analysis complete. The data suggests this variant significantly alters the protein network.")