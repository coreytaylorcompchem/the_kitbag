import requests
import pandas as pd
import time
import random

from tqdm import tqdm
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

PC_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def safe_get(url, max_retries=3, sleep_base=2):
    for attempt in range(max_retries):
        try:
            r = requests.get(url, timeout=20)
            if r.status_code == 200:
                return r.json()
        except Exception as e:
            logger.warning(f"Retry {attempt+1}/{max_retries} failed: {e}")
        time.sleep(sleep_base ** attempt + random.random())
    logger.error(f"❌ Failed GET after {max_retries} retries: {url}")
    return None


# ---------------------- NEW HELPERS ----------------------

def uniprot_to_gene_symbol(uniprot):
    """Resolve gene symbol from UniProt accession using UniProt REST API."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot}.json"
    js = safe_get(url)
    if js:
        genes = js.get("genes", [])
        for g in genes:
            if "geneName" in g:
                return g["geneName"]["value"]
    return None


def fetch_pubchem_aids(gene_id, gene_symbol):
    """
    Try all valid PubChem target namespaces that return AIDs:
    - geneid
    - genesymbol
    - proteinname
    """
    aids = set()
    endpoints = [
        f"{PC_BASE}/assay/target/geneid/{gene_id}/aids/JSON",
        f"{PC_BASE}/assay/target/genesymbol/{gene_symbol}/aids/JSON",
        f"{PC_BASE}/assay/target/proteinname/{gene_symbol}/aids/JSON",
    ]

    for url in endpoints:
        js = safe_get(url)
        if not js:
            continue
        idlist = js.get("IdentifierList", {})
        aids.update(idlist.get("AID", []))

    return sorted(aids)


def extract_pubchem_activity_rows(aid, assay_json, uniprot):
    """
    Extract full quantitative activity rows from the PC-AssayContainer
    data structure used by PubChem for real numeric readouts.
    """
    rows = []

    assay = (
        assay_json.get("Assay", {})
                  .get("PC_AssayContainer", [])
    )

    if not assay:
        return rows

    container = assay[0]  # 1 container per AID typically
    assay_data = container.get("assay", {})

    # 1) Lookup CIDs from sid → cid mapping
    sid_to_cid = {}
    for xref in assay_data.get("xref", []):
        sid = xref.get("sid")
        cid = xref.get("cid")
        if sid and cid:
            sid_to_cid[sid] = cid

    # 2) Find quantitative result columns
    result_cols = assay_data.get("results", [])

    for col in result_cols:
        heading = col.get("name")
        unit = col.get("unit")
        rows_data = col.get("data", [])   # each entry has { "sid": ..., "value": ... }

        # Skip non-numeric columns
        if not rows_data or "value" not in str(rows_data):
            continue

        for entry in rows_data:
            sid = entry.get("sid")
            value = entry.get("value")

            if sid not in sid_to_cid:
                continue

            cid = sid_to_cid[sid]

            rows.append({
                "AID": aid,
                "CID": cid,
                "ActivityType": heading,     # e.g., IC50, Ki, %Inhibition
                "ActivityValue": value,
                "ActivityUnit": unit,
                "uniprot_id": uniprot,
            })

    return rows


# ---------------------- MAIN TASK: Retrieve ----------------------

@register_task("retrieve_pubchem_bioactivities",
               category="Bioactivity",
               description="Retrieve bioactivity data from PubChem for a target.")
def retrieve_pubchem_bioactivities(config, data=None):

    uniprot = config.get("uniprot_id")
    readouts = config.get("readout", ["IC50"])
    if isinstance(readouts, str):
        readouts = [readouts]

    logger.info(f"🔍 Looking up GeneID for UniProt {uniprot}")

    # 1) UniProt → GeneID via NCBI Entrez
    gene_url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
        f"db=gene&term={uniprot}[Protein%20Accession]&retmode=json"
    )
    gene_json = safe_get(gene_url)

    gene_ids = []
    if gene_json and "esearchresult" in gene_json:
        gene_ids = gene_json["esearchresult"].get("idlist", [])

    # 2) If no GeneID, try to resolve gene symbol from config or UniProt
    gene_symbol = config.get("gene_symbol")
    if not gene_ids and not gene_symbol:
        gene_symbol = uniprot_to_gene_symbol(uniprot)
        if gene_symbol:
            logger.info(f"Resolved gene symbol '{gene_symbol}' from UniProt {uniprot}")
        else:
            logger.warning(f"No GeneID found for UniProt {uniprot} and no gene_symbol could be resolved.")
            return pd.DataFrame()

    # 3) If we have gene symbol but no GeneID, search Entrez by gene symbol
    if not gene_ids and gene_symbol:
        gene_symbol_url = (
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
            f"db=gene&term={gene_symbol}[Gene%20Name]&retmode=json"
        )
        gene_json2 = safe_get(gene_symbol_url)
        if gene_json2 and "esearchresult" in gene_json2:
            gene_ids = gene_json2["esearchresult"].get("idlist", [])

    if not gene_ids:
        logger.error(f"❌ Unable to find GeneID for {uniprot} / gene symbol {gene_symbol}")
        return pd.DataFrame()

    gene_id = gene_ids[0]
    logger.info(f"Found GeneID {gene_id} for UniProt {uniprot}")

    # 4) Retrieve PubChem assays
    aids = fetch_pubchem_aids(gene_id, gene_symbol)
    if not aids:
        logger.info(f"No PubChem assays found for GeneID {gene_id}")
        return pd.DataFrame()

    logger.info(f"Found {len(aids)} PubChem assays")

    rows = []
    for aid in tqdm(aids, desc="Fetching PubChem assays"):
        url = f"{PC_BASE}/assay/aid/{aid}/JSON"
        js = safe_get(url)
        if not js:
            continue
        rows.extend(extract_pubchem_activity_rows(aid, js, uniprot))

    if not rows:
        logger.info("No activity data found in PubChem assays.")
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    df = df.dropna(subset=["CID", "ActivityValue"])

    # Filter by readouts
    df["ActivityType_upper"] = df["ActivityType"].str.upper()
    df = df[df["ActivityType_upper"].isin([r.upper() for r in readouts])]
    df = df.drop(columns=["ActivityType_upper"])

    logger.info(f"Retrieved {len(df)} PubChem activity rows after filtering.")
    return df


# ---------------------- CLEAN ----------------------

@register_task("clean_pubchem_bioactivities",
               category="Bioactivity",
               description="Clean and standardize PubChem bioactivity records.")
def clean_pubchem_bioactivities(config, data):

    if data is None or not isinstance(data, pd.DataFrame):
        return {"df": pd.DataFrame(), "readout": None}

    df = data.copy()
    if df.empty:
        return {"df": df, "readout": None}

    # Standardize names
    df = df.rename(columns={
        "CID": "pubchem_cid",
        "ActivityType": "standard_type",
        "ActivityValue": "standard_value",
        "ActivityUnit": "standard_units",
    })

    # Force numeric activity values
    df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
    df = df.dropna(subset=["standard_value"])

    # Default units to nM
    df["standard_units"] = df["standard_units"].fillna("nM")

    readout = config.get("readout")
    return {"df": df, "readout": readout}


# ---------------------- SMILES ----------------------

@register_task("retrieve_pubchem_smiles",
               category="Bioactivity",
               description="Retrieve SMILES for PubChem CIDs.")
def retrieve_pubchem_smiles(config, data):

    if not isinstance(data, dict) or "df" not in data:
        logger.info("No DataFrame passed into retrieve_pubchem_smiles.")
        return {"df": pd.DataFrame(), "readout": None}

    df = data["df"]

    if df is None or not isinstance(df, pd.DataFrame) or df.empty:
        logger.info("Empty DataFrame → skipping SMILES retrieval.")
        return {"df": pd.DataFrame(), "readout": None}

    if "pubchem_cid" not in df.columns:
        logger.info("No 'pubchem_cid' column → skipping SMILES retrieval.")
        return {"df": df, "readout": data.get("readout")}

    cids = df["pubchem_cid"].dropna().unique()
    cids = [int(c) for c in cids]

    smiles_map = {}
    for cid in tqdm(cids, desc="Fetching PubChem SMILES"):
        url = f"{PC_BASE}/compound/cid/{cid}/property/CanonicalSMILES/JSON"
        js = safe_get(url)
        if js and "PropertyTable" in js:
            smiles = js["PropertyTable"]["Properties"][0].get("CanonicalSMILES")
            smiles_map[cid] = smiles

    df["smiles"] = df["pubchem_cid"].apply(lambda c: smiles_map.get(int(c), None))
    return {"df": df, "readout": data.get("readout")}
