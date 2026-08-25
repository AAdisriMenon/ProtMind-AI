import re
import requests
import streamlit as st
from Bio import Align

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

# ========= STEP 1: PARSING & VALIDATION =========
def parse_fasta(fasta_string: str) -> str:
    lines = fasta_string.strip().splitlines()
    sequence_lines = [line.strip() for line in lines if line and not line.startswith(">")]
    return "".join(sequence_lines)

def validate_sequence(seq: str) -> bool:
    clean_seq = "".join(seq.split())
    pattern = re.compile(r'^[ACDEFGHIKLMNPQRSTVWY]+$', re.IGNORECASE)
    return bool(pattern.match(clean_seq))

def validate_uniprot_id(uniprot_id: str) -> bool:
    pattern = re.compile(r'^[O,P,Q][0-9][A-Z,0-9]{3}[0-9]|[A-N,R-Z][0-9]([A-Z][A-Z,0-9]{2}[0-9]){1,2}$', re.IGNORECASE)
    return bool(pattern.match(uniprot_id.strip()))

def validate_mutation(mutation: str) -> bool:
    pattern = re.compile(r'^[ACDEFGHIKLMNPQRSTVWY]\d+[ACDEFGHIKLMNPQRSTVWY]$', re.IGNORECASE)
    return bool(pattern.match(mutation.strip()))

# ========= STEP 2: DATA RETRIEVAL =========
def fetch_uniprot_data(uniprot_id: str) -> dict:
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
    response = requests.get(url)
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
            "Sequence Length": data.get("sequence", {}).get("length", 0)
        }
    return {"Error": f"Could not retrieve data for {uniprot_id} (Status: {response.status_code})"}

def fetch_alphafold_url(uniprot_id: str) -> str:
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    response = requests.get(url)
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
    is_valid_pos = (position <= len(sequence)) and (sequence[position-1].upper() == wt_aa)

    return {
        "Mutation": mutation,
        "Wildtype AA": wt_aa,
        "Mutant AA": mut_aa,
        "Position": position,
        "Positional Match": is_valid_pos,
        "WT Property": get_amino_acid_property(wt_aa),
        "Mutant Property": get_amino_acid_property(mut_aa),
        "Property Shift": f"{get_amino_acid_property(wt_aa)} ➡️ {get_amino_acid_property(mut_aa)}"
    }

def calculate_alignment_score(wildtype_seq: str, mutant_seq: str) -> float:
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    alignments = aligner.align(wildtype_seq, mutant_seq)
    best_alignment = alignments[0]
    return (best_alignment.score / len(wildtype_seq)) * 100

# ========= STEP 4: MUTATION ANALYSIS & FUNCTIONAL PROTEOMICS =========

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
            
            # Extract the top interacting proteins based on the limit
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

def predict_pathogenicity(mutation: str) -> dict:
    """
    A framework to aggregate pathogenicity scores. 
    (In a production environment, this would query the Ensembl VEP or UniProt Variation API).
    """
    # For now, we simulate the classifier engine based on your workflow diagram
    return {
        "SIFT": {"score": 0.02, "prediction": "Deleterious"},
        "PolyPhen-2": {"score": 0.95, "prediction": "Probably Damaging"},
        "AlphaMissense": {"score": 0.88, "prediction": "Pathogenic"}
    }

# ========= STEP 5: STRUCTURAL PROTEOMICS =========

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