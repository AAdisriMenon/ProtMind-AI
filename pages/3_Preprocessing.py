import json
import streamlit as st
import pandas as pd
from backend import (
    apply_custom_css,
    render_progress,
    calculate_alignment_score,
    extract_features,
    extract_sequence_profile,
    fetch_pfam_domains,
    scan_sequence_motifs,
)

st.set_page_config(
    page_title="Step 03 | ProtMind AI", page_icon="⚙️", layout="wide"
)
apply_custom_css()
render_progress(3)

st.markdown(
    "<h2>⚙️ Step 03: Preprocessing & Motif/Domain Search</h2>",
    unsafe_allow_html=True,
)

# Guard clause for payload
if "payload" not in st.session_state or not st.session_state["payload"]:
    st.warning(
        "⚠️ No active protein data found. Please complete Step 01: User Input & Validation first."
    )
    st.page_link("pages/1_User_Input.py", label="⬅️ Back to Step 01", icon="⚡")
    st.stop()

payload = st.session_state["payload"]
uniprot_id = payload.get("uniprot_id", "").strip()
raw_seq = payload.get("sequence", "").strip()
mutation = payload.get("mutation", "").strip()

if not raw_seq:
    st.error(
        "❌ Sequence is missing. Please ensure sequence retrieval or entry succeeded in previous steps."
    )
    st.page_link("pages/2_Data_Retrieval.py", label="⬅️ Back to Step 02", icon="🔍")
    st.stop()

retrieved = st.session_state.get("retrieved_data", {})
uniprot_metadata = retrieved.get("metadata", {})
uniprot_features = uniprot_metadata.get("Features", []) if isinstance(uniprot_metadata, dict) else []

cache_key = f"{uniprot_id}:{mutation}:{len(raw_seq)}"

# Cache preprocessing, Pfam domain search, and motif scan in session_state
if "preprocessed_data" not in st.session_state or st.session_state.get("processed_key") != cache_key:
    with st.spinner("⏳ Extracting sequence features, scanning motifs, and querying Pfam domains..."):
        try:
            features = extract_features(raw_seq, mutation)
        except ValueError:
            st.error(f"❌ Could not parse mutation position from '{mutation}'. Expected format like R75H or R175H.")
            st.stop()

        if features.get("Positional Match"):
            mutant_list = list(raw_seq)
            mutant_list[features["Position"] - 1] = features["Mutant AA"]
            mutant_sequence = "".join(mutant_list)
        else:
            mutant_sequence = raw_seq

        alignment_score = calculate_alignment_score(raw_seq, mutant_sequence)
        seq_profile = extract_sequence_profile(raw_seq)
        pfam_domains = fetch_pfam_domains(
            uniprot_id,
            mutation_pos=features.get("Position"),
            uniprot_metadata=uniprot_metadata if isinstance(uniprot_metadata, dict) else None,
        )
        motif_results = scan_sequence_motifs(
            raw_seq,
            mutation_pos=features.get("Position"),
            uniprot_features=uniprot_features,
        )

        st.session_state["preprocessed_data"] = {
            "features": features,
            "mutant_sequence": mutant_sequence,
            "alignment_score": alignment_score,
            "seq_profile": seq_profile,
            "pfam_domains": pfam_domains,
            "motif_results": motif_results,
        }
        st.session_state["processed_key"] = cache_key

proc_data = st.session_state["preprocessed_data"]
features = proc_data["features"]
alignment_score = proc_data["alignment_score"]
seq_profile = proc_data["seq_profile"]
pfam_domains = proc_data["pfam_domains"]
motif_results = proc_data["motif_results"]

# ================= TABS FOR STEP 03 =================
tab_seq, tab_pfam, tab_motif = st.tabs([
    "🧪 Sequence Features & Mutation Check",
    "🧬 Conserved Domains (Pfam)",
    "🔎 Motif & Functional Site Search",
])

with tab_seq:
    col_feat, col_align = st.columns(2)

    with col_feat:
        st.markdown("<h3>Mutation & Physicochemical Shift</h3>", unsafe_allow_html=True)

        if not features.get("In Range"):
            st.error(
                f"⚠️ Mutation position {features.get('Position')} is outside the sequence "
                f"(length {len(raw_seq)}). Double check the mutation and sequence match."
            )
        elif not features.get("Positional Match"):
            st.error(
                f"⚠️ Positional mismatch: Expected `{features.get('Wildtype AA')}` at position {features.get('Position')}, "
                f"but found `{raw_seq[features['Position'] - 1]}`."
            )
        else:
            st.success("✅ Mutation position verified against wildtype sequence.")

        f_col1, f_col2 = st.columns(2)
        with f_col1:
            st.metric("Wildtype Residue", f"{features['Wildtype AA']} ({features['WT Property']})")
        with f_col2:
            st.metric("Mutant Residue", f"{features['Mutant AA']} ({features['Mutant Property']})")

        if features["WT Property"] != features["Mutant Property"]:
            st.warning(f"**Biophysical Shift:** {features['Property Shift']}")
        else:
            st.info(f"**Conserved Class:** Maintained as {features['WT Property']}")

    with col_align:
        st.markdown("<h3>Sequence Alignment & Global Profile</h3>", unsafe_allow_html=True)
        st.metric(label="WT vs Mutant Identity Score", value=f"{alignment_score:.2f}%")
        st.progress(min(max(alignment_score / 100.0, 0.0), 1.0))

        m1, m2, m3 = st.columns(3)
        m1.metric("Molecular Weight", f"{seq_profile.get('Molecular Weight (Da)', 0):,.1f} Da")
        m2.metric("Theoretical pI", f"{seq_profile.get('Isoelectric Point (pI)', 0)}")
        m3.metric("GRAVY", f"{seq_profile.get('GRAVY (Hydropathicity)', 0)}")

        m4, m5, m6 = st.columns(3)
        m4.metric("Aromaticity", f"{seq_profile.get('Aromaticity', 0)}")
        m5.metric("Instability Index", f"{seq_profile.get('Instability Index', 0)}")
        m6.metric("Charge (pH 7.0)", f"{seq_profile.get('Charge at pH 7.0', 0)}")

    st.markdown("#### 📊 Amino Acid Class Composition (%)")
    comp = seq_profile.get("Class Composition (%)", {})
    if comp:
        st.bar_chart(comp)

    with st.expander("📄 View Extracted Feature Dictionaries"):
        st.json({"mutation_features": features, "sequence_profile": seq_profile})

with tab_pfam:
    st.markdown("<h3>Conserved Domains (Pfam)</h3>", unsafe_allow_html=True)
    if pfam_domains:
        df_pfam = pd.DataFrame(pfam_domains)
        st.dataframe(df_pfam, use_container_width=True, hide_index=True)

        hit_domains = [d for d in pfam_domains if d.get("Mutation in Domain") == "Yes"]
        if hit_domains:
            names = ", ".join(f"{d['Pfam ID']} ({d['Domain Name']})" for d in hit_domains)
            st.warning(f"🎯 Target mutation **{mutation}** falls inside conserved Pfam domain(s): **{names}**")
        else:
            st.info(f"ℹ️ Target mutation **{mutation}** does not fall within the mapped Pfam domain boundaries.")
    else:
        st.warning("⚠️ No conserved Pfam domains were returned for this accession.")

with tab_motif:
    st.markdown("<h3>Motifs & Functional Sites</h3>", unsafe_allow_html=True)

    curated = motif_results.get("curated_motifs", [])
    pattern_hits = motif_results.get("pattern_hits", [])

    st.markdown("#### 1️⃣ UniProtKB Curated Motifs & Functional Sites")
    if curated:
        st.dataframe(pd.DataFrame(curated), use_container_width=True, hide_index=True)
        curated_overlaps = [m for m in curated if m.get("Overlaps Mutation") == "Yes"]
        if curated_overlaps:
            st.warning(
                f"🎯 Mutation **{mutation}** directly overlaps curated site(s): "
                + ", ".join(f"{m['Type']} ({m['Description']})" for m in curated_overlaps)
            )
    else:
        st.info("No curated Motif / Active site / Binding site annotations found in UniProtKB for this entry.")

    st.markdown("#### 2️⃣ Sequence Pattern Motif Scan (PROSITE-style)")
    if pattern_hits:
        st.dataframe(pd.DataFrame(pattern_hits), use_container_width=True, hide_index=True)
        pat_overlaps = [p for p in pattern_hits if p.get("Overlaps Mutation") == "Yes"]
        if pat_overlaps:
            st.warning(
                f"🎯 Mutation **{mutation}** overlaps predicted sequence motif(s): "
                + ", ".join(f"{p['Pattern ID']} ({p['Motif']} @ {p['Start']}-{p['End']})" for p in pat_overlaps)
            )
    else:
        st.info("No PROSITE-style sequence motif patterns matched this sequence.")

st.markdown("---")
st.success("🎉 **3-Step Prototype Complete!** You have completed User Input & Validation, Data Retrieval, and Preprocessing & Motif/Domain Search.")

export_summary = {
    "step_01_input": {
        "uniprot_id": uniprot_id,
        "mutation": mutation,
        "source": payload.get("source"),
        "sequence_length": len(raw_seq),
    },
    "step_02_retrieval": {
        "protein_name": uniprot_metadata.get("Protein Name") if isinstance(uniprot_metadata, dict) else None,
        "gene": uniprot_metadata.get("Gene") if isinstance(uniprot_metadata, dict) else None,
        "organism": uniprot_metadata.get("Organism") if isinstance(uniprot_metadata, dict) else None,
        "alphafold_pdb_url": retrieved.get("pdb_url"),
    },
    "step_03_preprocessing": {
        "mutation_features": features,
        "sequence_profile": seq_profile,
        "pfam_domains": pfam_domains,
        "motif_results": motif_results,
    },
}

col_dl, col_home = st.columns(2)
with col_dl:
    st.download_button(
        label="📥 Download 3-Step Analysis Summary (JSON)",
        data=json.dumps(export_summary, indent=2),
        file_name=f"{uniprot_id}_{mutation}_protmind_3step.json",
        mime="application/json",
    )
with col_home:
    st.page_link("app.py", label="🔄 Return to Home / Start New Analysis", icon="🧠")
