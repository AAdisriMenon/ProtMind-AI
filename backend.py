import re
import requests
import streamlit as st
from Bio import Align
from Bio.Align import substitution_matrices

# ========= UI / CSS STYLING =========
def apply_custom_css():
    st.markdown("""
    <style>
    .stApp { background: linear-gradient(135deg, black 100%, black 100%, black 100%); background-attachment: fixed; }
    .block-container { background: rgba(22, 26, 35, 0.75); backdrop-filter: blur(10px); border: 1px solid rgba(0, 245, 212, 0.25); border-radius: 40px; padding: 2.5rem 3rem!important; margin-top: 2rem; box-shadow: 0 8px 32px 0 rgba(0, 0, 0, 0.6); }
    h1 { background: linear-gradient(90deg, #00F5D4, #9B5DE5); -webkit-background-clip: text; -webkit-text-fill-color: transparent; font-weight: 900!important; text-align: center; font-size: 3.5rem!important; margin-bottom: 1.5rem!important; }
    h2, h3 { color: #FFFFFF!important; text-align: center; font-weight: 800!important; }
    p, label, div, span { color: #F8F9FA!important; font-weight: 500!important; }
    label { color: #FFFFFF!important; font-weight: 700!important; }
    .stTextInput>div>div>input,.stTextArea textarea { background-color: rgba(15, 30, 60, 0.9)!important; border: 1.5px solid #00F5D4!important; color: #FFFFFF!important; border-radius: 12px!important; font-weight: 600!important; }
    .stTextInput>div>div>input:focus,.stTextArea textarea:focus { border: 2px solid #9B5DE5!important; box-shadow: 0 0 0 2px #9B5DE5!important; }
    .stButton>button { background: linear-gradient(90deg, #00F5D4 0%, #9B5DE5 100%); color: #8B5CF6!important; font-weight: 800; border: none; border-radius: 14px; padding: 0.8rem 2rem; width: 100%; transition: all 0.3s ease; margin-top: 1rem; }
    .stButton>button:hover { transform: translateY(-3px); box-shadow: 0 0 30px rgba(0, 245, 212, 0.8); }
    [data-testid="stHeader"] { background-color: rgba(0,0,0,0); }
    hr { border-color: rgba(0, 245, 212, 0.3)!important; }
    .stJson { background-color: rgba(10, 25, 47, 0.9)!important; border-radius: 12px!important; border: 1px solid #00F5D4!important; }
    [data-testid="stFileUploadDropzone"] { background-color: rgba(15, 30, 60, 0.9)!important; border: 1.5px dashed #00F5D4!important; border-radius: 12px!important; }
    </style>
    """, unsafe_allow_html=True)


# ========= PIPELINE PROGRESS TRACKER =========
PIPELINE_STEPS = [
    "User Input",
    "Data Retrieval",
    "Preprocessing",
    "Mutation Analysis",
    "Structure Prediction",
]


def render_progress(current_step: int):
    """Renders a sidebar checklist showing where the user is in the 5-step pipeline."""
    st.sidebar.markdown("### 🧭 Pipeline Progress")
    for i, name in enumerate(PIPELINE_STEPS, start=1):
        if i < current_step:
            st.sidebar.markdown(f"✅ {i}. {name}")
        elif i == current_step:
            st.sidebar.markdown(f"➡️ **{i}. {name}**")
        else:
            st.sidebar.markdown(f"⬜ {i}. {name}")
    st.sidebar.divider()
    if st.sidebar.button("🔄 Reset Session"):
        for key in list(st.session_state.keys()):
            del st.session_state[key]
        st.switch_page("app.py")


# ========= STEP 1: PARSING & VALIDATION =========
def parse_fasta(fasta_string: str) -> str:
    lines = fasta_string.strip().splitlines()
    sequence_lines = [line.strip() for line in lines if line and not line.startswith(">")]
    return "".join(sequence_lines)


def validate_sequence(seq: str) -> bool:
    clean_seq = "".join(seq.split())
    if not clean_seq:
        return False
    pattern = re.compile(r'^[ACDEFGHIKLMNPQRSTVWY]+$', re.IGNORECASE)
    return bool(pattern.match(clean_seq))


def validate_uniprot_id(uniprot_id: str) -> bool:
    # Fixed: original regex leaked the alternation outside the anchors
    # (^ only bound the first alternative, $ only bound the second),
    # so it silently accepted malformed IDs. Now both branches are anchored.
    pattern = re.compile(
        r'^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$',
        re.IGNORECASE,
    )
    return bool(pattern.match(uniprot_id.strip()))


def validate_mutation(mutation: str) -> bool:
    pattern = re.compile(r'^[ACDEFGHIKLMNPQRSTVWY]\d+[ACDEFGHIKLMNPQRSTVWY]$', re.IGNORECASE)
    return bool(pattern.match(mutation.strip()))


# ========= STEP 2: DATA RETRIEVAL =========
@st.cache_data(ttl=3600, show_spinner=False)
def fetch_uniprot_data(uniprot_id: str) -> dict:
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
    try:
        response = requests.get(url, timeout=10)
    except requests.exceptions.RequestException as exc:
        return {"Error": f"Network error contacting UniProt: {exc}"}

    if response.status_code == 200:
        data = response.json()
        try:
            protein_name = data.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "Unknown")
        except AttributeError:
            protein_name = "Unknown"
        return {
            "Entry": data.get("primaryAccession", uniprot_id),
            "Protein Name": protein_name,
            "Gene": data.get("genes", [{}])[0].get("geneName", {}).get("value", "Unknown"),
            "Organism": data.get("organism", {}).get("scientificName", "Unknown"),
            "Sequence Length": data.get("sequence", {}).get("length", 0),
            # Fixed: the live sequence text was never actually retrieved, so the
            # "auto-fill sequence" logic in Step 2 was always a dead branch.
            "Sequence": data.get("sequence", {}).get("value", ""),
        }
    return {"Error": f"Could not retrieve data for {uniprot_id} (Status: {response.status_code})"}


@st.cache_data(ttl=3600, show_spinner=False)
def fetch_alphafold_url(uniprot_id: str) -> str:
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        response = requests.get(url, timeout=10)
    except requests.exceptions.RequestException:
        return None

    if response.status_code == 200:
        data = response.json()
        if len(data) > 0:
            return data[0].get("pdbUrl", None)
    return None


# ========= STEP 3: PREPROCESSING & BIOPHYSICS =========
def get_amino_acid_property(aa: str) -> str:
    hydrophobic, polar, positive, negative = {'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'}, {'S', 'T', 'N', 'Q', 'C'}, {'R', 'H', 'K'}, {'D', 'E'}
    if aa in hydrophobic: return "Hydrophobic"
    if aa in polar: return "Polar/Neutral"
    if aa in positive: return "Positively Charged"
    if aa in negative: return "Negatively Charged"
    return "Special/Other"


def extract_features(sequence: str, mutation: str) -> dict:
    wt_aa = mutation[0].upper()
    mut_aa = mutation[-1].upper()
    position = int(mutation[1:-1])

    # Fixed: this used to index straight into the sequence, which threw an
    # uncaught IndexError (crashing the page) whenever the mutation position
    # was out of range. Now the bounds check happens before indexing.
    is_in_range = 1 <= position <= len(sequence)
    is_valid_pos = is_in_range and sequence[position - 1].upper() == wt_aa

    return {
        "Mutation": mutation,
        "Wildtype AA": wt_aa,
        "Mutant AA": mut_aa,
        "Position": position,
        "Positional Match": is_valid_pos,
        "In Range": is_in_range,
        "WT Property": get_amino_acid_property(wt_aa),
        "Mutant Property": get_amino_acid_property(mut_aa),
        "Property Shift": f"{get_amino_acid_property(wt_aa)} ➡️ {get_amino_acid_property(mut_aa)}"
    }


def calculate_alignment_score(wildtype_seq: str, mutant_seq: str) -> float:
    if not wildtype_seq:
        return 0.0
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    alignments = aligner.align(wildtype_seq, mutant_seq)
    best_alignment = alignments[0]
    return (best_alignment.score / len(wildtype_seq)) * 100


# ========= STEP 4: MUTATION ANALYSIS & FUNCTIONAL PROTEOMICS =========
@st.cache_data(ttl=3600, show_spinner=False)
def fetch_string_ppi(uniprot_id: str, limit: int = 5) -> list:
    """
    Fetches the top protein-protein interactions from the STRING database.
    Assumes Human species (Taxonomy ID: 9606).
    """
    url = f"https://string-db.org/api/json/network?identifiers={uniprot_id}&species=9606"

    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            interactions = []
            for i, edge in enumerate(data):
                if i >= limit:
                    break
                interactions.append({
                    "Target": edge.get("preferredName_B", "Unknown"),
                    "Interaction Score": edge.get("score", 0.0)
                })
            return interactions
    except requests.exceptions.RequestException:
        return []

    return []


_BLOSUM62 = substitution_matrices.load("BLOSUM62")


def predict_pathogenicity(mutation: str) -> dict:
    """
    Upgraded: this used to return three hardcoded numbers regardless of input.
    It now derives a real, explainable risk score from the BLOSUM62
    substitution matrix (how conservative/disruptive the swap is) plus a
    penalty if the amino acid's physicochemical class changes.

    This is still a fast heuristic, not a clinical-grade predictor -- a
    production build would call real classifiers (SIFT, PolyPhen-2,
    AlphaMissense, or the Ensembl VEP API) instead.
    """
    wt_aa, mut_aa = mutation[0].upper(), mutation[-1].upper()

    try:
        raw_score = _BLOSUM62[wt_aa, mut_aa]
    except KeyError:
        raw_score = 0

    # BLOSUM62 scores roughly range from -4 (very disruptive) to +11 (identical).
    # Map that onto a 0 (benign) -> 1 (deleterious) risk scale.
    normalized_risk = max(0.0, min(1.0, (4 - raw_score) / 8))

    wt_property = get_amino_acid_property(wt_aa)
    mut_property = get_amino_acid_property(mut_aa)
    if wt_property != mut_property:
        normalized_risk = min(1.0, normalized_risk + 0.15)

    sift_score = round(max(0.0, 1 - normalized_risk), 3)
    polyphen_score = round(normalized_risk, 3)
    alpha_missense_score = round(min(1.0, (normalized_risk * 0.7) + (polyphen_score * 0.3)), 3)

    def sift_label(s):
        return "Deleterious" if s < 0.05 else "Tolerated"

    def polyphen_label(s):
        if s > 0.85:
            return "Probably Damaging"
        if s > 0.45:
            return "Possibly Damaging"
        return "Benign"

    def alpha_label(s):
        if s > 0.7:
            return "Pathogenic"
        if s > 0.4:
            return "Ambiguous"
        return "Likely Benign"

    return {
        "SIFT": {"score": sift_score, "prediction": sift_label(sift_score)},
        "PolyPhen-2": {"score": polyphen_score, "prediction": polyphen_label(polyphen_score)},
        "AlphaMissense": {"score": alpha_missense_score, "prediction": alpha_label(alpha_missense_score)},
    }


# ========= STEP 5: STRUCTURAL PROTEOMICS =========
@st.cache_data(ttl=3600, show_spinner=False)
def fetch_pdb_content(pdb_url: str) -> str:
    """Fetches the raw PDB file text from the AlphaFold database URL."""
    if not pdb_url:
        return ""
    try:
        response = requests.get(pdb_url, timeout=10)
        if response.status_code == 200:
            return response.text
    except requests.exceptions.RequestException:
        return ""
    return ""
