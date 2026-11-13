import re
import requests
import streamlit as st
from Bio import Align
from Bio.SeqUtils.ProtParam import ProteinAnalysis

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


# ========= PIPELINE PROGRESS TRACKER (3-STEP PROTOTYPE) =========
PIPELINE_STEPS = [
    "User Input & Validation",
    "Data Retrieval",
    "Preprocessing & Motif/Domain Search",
]


def render_progress(current_step: int):
    """Renders a sidebar checklist showing where the user is in the 3-step prototype pipeline."""
    st.sidebar.markdown("### 🧭 Pipeline Progress")
    for i, name in enumerate(PIPELINE_STEPS, start=1):
        if i < current_step:
            st.sidebar.markdown(f"✅ Step 0{i}: {name}")
        elif i == current_step:
            st.sidebar.markdown(f"➡️ **Step 0{i}: {name}**")
        else:
            st.sidebar.markdown(f"⬜ Step 0{i}: {name}")
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
    """Fetches protein annotations, sequence, features, and Pfam cross-refs from UniProtKB."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
    try:
        response = requests.get(url, timeout=10)
    except requests.exceptions.RequestException as exc:
        return {"Error": f"Network error contacting UniProt: {exc}"}

    if response.status_code != 200:
        return {"Error": f"Could not retrieve data for {uniprot_id} (Status: {response.status_code})"}

    data = response.json()
    try:
        protein_name = (
            data.get("proteinDescription", {})
            .get("recommendedName", {})
            .get("fullName", {})
            .get("value", "Unknown")
        )
    except AttributeError:
        protein_name = "Unknown"

    # Extract function and subcellular location comments
    function_texts = []
    subcellular_locations = []
    for comment in data.get("comments", []):
        ctype = comment.get("commentType", "")
        if ctype == "FUNCTION":
            for t in comment.get("texts", []):
                if t.get("value"):
                    function_texts.append(t["value"])
        elif ctype == "SUBCELLULAR LOCATION":
            for loc in comment.get("subcellularLocations", []):
                val = loc.get("location", {}).get("value")
                if val:
                    subcellular_locations.append(val)

    # Extract annotated domains, motifs, and functional sites
    relevant_types = {"Domain", "Region", "Motif", "Active site", "Binding site", "Zinc finger", "DNA binding", "Site"}
    annotated_features = []
    for feat in data.get("features", []):
        ftype = feat.get("type", "")
        if ftype in relevant_types:
            loc = feat.get("location", {})
            start = loc.get("start", {}).get("value")
            end = loc.get("end", {}).get("value")
            annotated_features.append({
                "Type": ftype,
                "Description": feat.get("description", "—") or "—",
                "Start": start,
                "End": end,
            })

    # Extract Pfam cross-references
    pfam_refs = []
    for xref in data.get("uniProtKBCrossReferences", []):
        if xref.get("database") == "Pfam":
            props = {p.get("key"): p.get("value") for p in xref.get("properties", [])}
            pfam_refs.append({
                "Pfam ID": xref.get("id", ""),
                "Domain Name": props.get("EntryName", "Unknown"),
            })

    keywords = [kw.get("name") for kw in data.get("keywords", []) if kw.get("name")]

    return {
        "Entry": data.get("primaryAccession", uniprot_id),
        "Protein Name": protein_name,
        "Gene": data.get("genes", [{}])[0].get("geneName", {}).get("value", "Unknown") if data.get("genes") else "Unknown",
        "Organism": data.get("organism", {}).get("scientificName", "Unknown"),
        "Sequence Length": data.get("sequence", {}).get("length", 0),
        "Sequence": data.get("sequence", {}).get("value", ""),
        "Function": " ".join(function_texts) if function_texts else "No functional summary available.",
        "Subcellular Location": ", ".join(dict.fromkeys(subcellular_locations)) if subcellular_locations else "Not specified",
        "Keywords": keywords[:12],
        "Features": annotated_features,
        "Pfam Cross-Refs": pfam_refs,
    }


@st.cache_data(ttl=3600, show_spinner=False)
def fetch_alphafold_prediction(uniprot_id: str) -> dict:
    """Fetches AlphaFold DB prediction metadata including PDB/CIF URLs, residue coverage, and fragment info."""
    clean_accession = uniprot_id.split("-")[0].strip()
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{clean_accession}"
    try:
        response = requests.get(api_url, timeout=10)
    except requests.exceptions.RequestException:
        return {}

    if response.status_code == 200:
        data = response.json()
        if not data or not isinstance(data, list):
            return {}

        chosen = data[0]
        return {
            "pdbUrl": chosen.get("pdbUrl"),
            "cifUrl": chosen.get("cifUrl"),
            "entryId": chosen.get("entryId"),
            "uniprotStart": chosen.get("uniprotStart", 1),
            "uniprotEnd": chosen.get("uniprotEnd"),
            "uniprotSequence": chosen.get("uniprotSequence"),
            "latestVersion": chosen.get("latestVersion"),
            "allFragments": data,
        }
    return {}


@st.cache_data(ttl=3600, show_spinner=False)
def fetch_alphafold_url(uniprot_id: str) -> str:
    pred = fetch_alphafold_prediction(uniprot_id)
    return pred.get("pdbUrl")


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


# ========= STEP 3: PREPROCESSING & MOTIF/DOMAIN SEARCH =========
def get_amino_acid_property(aa: str) -> str:
    hydrophobic = {'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'}
    polar = {'S', 'T', 'N', 'Q', 'C'}
    positive = {'R', 'H', 'K'}
    negative = {'D', 'E'}
    if aa in hydrophobic:
        return "Hydrophobic"
    if aa in polar:
        return "Polar/Neutral"
    if aa in positive:
        return "Positively Charged"
    if aa in negative:
        return "Negatively Charged"
    return "Special/Other"


def extract_features(sequence: str, mutation: str) -> dict:
    wt_aa = mutation[0].upper()
    mut_aa = mutation[-1].upper()
    position = int(mutation[1:-1])

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
        "Property Shift": f"{get_amino_acid_property(wt_aa)} ➡️ {get_amino_acid_property(mut_aa)}",
    }


def extract_sequence_profile(sequence: str) -> dict:
    """Computes global biophysical sequence features using Biopython ProteinAnalysis."""
    clean_seq = "".join(sequence.split()).upper()
    if not clean_seq:
        return {}

    analyser = ProteinAnalysis(clean_seq)
    total = len(clean_seq)
    counts = {
        "Hydrophobic": sum(1 for aa in clean_seq if aa in {'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'}),
        "Polar/Neutral": sum(1 for aa in clean_seq if aa in {'S', 'T', 'N', 'Q', 'C'}),
        "Positively Charged": sum(1 for aa in clean_seq if aa in {'R', 'H', 'K'}),
        "Negatively Charged": sum(1 for aa in clean_seq if aa in {'D', 'E'}),
        "Special (G, P)": sum(1 for aa in clean_seq if aa in {'G', 'P'}),
    }
    composition_pct = {k: round((v / total) * 100, 2) for k, v in counts.items()}

    return {
        "Length (aa)": total,
        "Molecular Weight (Da)": round(analyser.molecular_weight(), 2),
        "Isoelectric Point (pI)": round(analyser.isoelectric_point(), 2),
        "GRAVY (Hydropathicity)": round(analyser.gravy(), 3),
        "Aromaticity": round(analyser.aromaticity(), 3),
        "Instability Index": round(analyser.instability_index(), 2),
        "Charge at pH 7.0": round(analyser.charge_at_pH(7.0), 2),
        "Class Composition (%)": composition_pct,
    }


def calculate_alignment_score(wildtype_seq: str, mutant_seq: str) -> float:
    if not wildtype_seq:
        return 0.0
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    alignments = aligner.align(wildtype_seq, mutant_seq)
    best_alignment = alignments[0]
    return (best_alignment.score / len(wildtype_seq)) * 100


@st.cache_data(ttl=3600, show_spinner=False)
def fetch_pfam_domains(uniprot_id: str, mutation_pos: int = None, uniprot_metadata: dict = None) -> list:
    """
    Queries the InterPro/Pfam REST API for conserved Pfam domains and their residue coordinates.
    Falls back to UniProtKB Pfam cross-references if the InterPro API is unavailable.
    """
    clean_id = uniprot_id.split("-")[0].strip()
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{clean_id}/?page_size=50"
    domains = []
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            results = response.json().get("results", [])
            for item in results:
                meta = item.get("metadata", {})
                pfam_id = meta.get("accession", "")
                name = meta.get("name", "Unknown Domain")
                for prot in item.get("proteins", []):
                    for loc in prot.get("entry_protein_locations", []):
                        for frag in loc.get("fragments", []):
                            start = frag.get("start")
                            end = frag.get("end")
                            contains_mut = (
                                "Yes" if (mutation_pos and start and end and start <= mutation_pos <= end)
                                else ("No" if mutation_pos else "—")
                            )
                            domains.append({
                                "Pfam ID": pfam_id,
                                "Domain Name": name,
                                "Start": start,
                                "End": end,
                                "Length": (end - start + 1) if (start and end) else None,
                                "Mutation in Domain": contains_mut,
                                "Source": "InterPro / Pfam API",
                            })
    except requests.exceptions.RequestException:
        pass

    # Fallback to UniProtKB Pfam cross-references if InterPro returned nothing
    if not domains and uniprot_metadata:
        for xref in uniprot_metadata.get("Pfam Cross-Refs", []):
            domains.append({
                "Pfam ID": xref.get("Pfam ID", ""),
                "Domain Name": xref.get("Domain Name", "Unknown"),
                "Start": None,
                "End": None,
                "Length": None,
                "Mutation in Domain": "—",
                "Source": "UniProtKB Cross-Ref",
            })

    return domains


# Common PROSITE-style sequence motif patterns
_SEQUENCE_MOTIF_PATTERNS = [
    ("PS00001", "N-glycosylation site", r"N[^P][ST][^P]"),
    ("PS00004", "cAMP/cGMP-dependent protein kinase phosphorylation site", r"[RK]{2}.[ST]"),
    ("PS00005", "Protein kinase C (PKC) phosphorylation site", r"[ST].[RK]"),
    ("PS00006", "Casein kinase II (CK2) phosphorylation site", r"[ST].{2}[DE]"),
    ("PS00007", "Tyrosine kinase phosphorylation site", r"[RK].{2,3}[DE].{2,3}Y"),
    ("PS00008", "N-myristoylation site", r"G[^EDRKHPFYW].{2}[STAGCN][^P]"),
    ("PS00016", "RGD cell attachment sequence", r"RGD"),
    ("PS00028", "Zinc finger C2H2-type signature", r"C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H"),
]


def scan_sequence_motifs(sequence: str, mutation_pos: int = None, uniprot_features: list = None) -> dict:
    """
    Identifies functional motifs and sites from:
    1) UniProtKB curated annotations (Motif, Active site, Binding site, Zinc finger, DNA binding, Site)
    2) PROSITE-style sequence pattern scanning
    """
    clean_seq = "".join(sequence.split()).upper()

    # 1. UniProt curated motifs & sites
    curated_motifs = []
    motif_types = {"Motif", "Active site", "Binding site", "Zinc finger", "DNA binding", "Site"}
    for feat in (uniprot_features or []):
        if feat.get("Type") in motif_types:
            start = feat.get("Start")
            end = feat.get("End")
            overlaps = (
                "Yes" if (mutation_pos and start and end and start <= mutation_pos <= end)
                else ("No" if mutation_pos else "—")
            )
            curated_motifs.append({
                "Type": feat.get("Type"),
                "Description": feat.get("Description", "—"),
                "Start": start,
                "End": end,
                "Overlaps Mutation": overlaps,
            })

    # 2. Sequence regex motif scan
    pattern_hits = []
    for prosite_id, motif_name, regex in _SEQUENCE_MOTIF_PATTERNS:
        for match in re.finditer(f"(?=({regex}))", clean_seq):
            start = match.start(1) + 1
            end = match.end(1)
            overlaps = (
                "Yes" if (mutation_pos and start <= mutation_pos <= end)
                else ("No" if mutation_pos else "—")
            )
            pattern_hits.append({
                "Pattern ID": prosite_id,
                "Motif": motif_name,
                "Start": start,
                "End": end,
                "Matched Sequence": match.group(1),
                "Overlaps Mutation": overlaps,
            })

    return {
        "curated_motifs": curated_motifs,
        "pattern_hits": pattern_hits,
    }
