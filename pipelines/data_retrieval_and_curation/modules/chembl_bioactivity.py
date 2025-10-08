from pipeline.task_registry import register_task
from chembl_webresource_client.new_client import new_client

import math
import time
import random

import pandas as pd
from tqdm import tqdm

from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

def retry_activity_query(filters, max_retries=3):
    activity = new_client.activity
    for attempt in range(max_retries):
        try:
            return activity.filter(**filters)
        except Exception as e:
            logger.warning(f"Retry {attempt+1}/{max_retries} failed: {e}")
            time.sleep(2 ** attempt + random.random())
    logger.error("❌ Final retry failed.")
    return []

@register_task("retrieve_chembl_bioactivities", 
               category='Bioactivity',
               description="Retrieve bioactivity data from CHEMBL.")
def retrieve_chembl_bioactivities(config, data=None):
    uniprot_id = config.get("uniprot_id")
    assay_type = config.get("assay_type")
    relation = config.get("relation")
    readouts = config.get("readout", ["IC50"])
    
    if isinstance(readouts, str):
        readouts = [readouts]

    target = new_client.target
    # activity = new_client.activity

    # Lookup target by uniprot
    target_query = target.filter(target_components__accession=uniprot_id)
    if not target_query:
        logger.info(f"❌ No targets found for UniProt ID: {uniprot_id}")
        return {"df": pd.DataFrame(), "readout": None}

    target_chembl_id = target_query[0]["target_chembl_id"]

    # Base filter
    filters = {
        "target_chembl_id": target_chembl_id,
        "standard_type__in": readouts
    }

    if relation:
        filters["standard_relation"] = relation

    if assay_type:
        filters["assay_type"] = assay_type

    # Query bioactivities
    activities = retry_activity_query(filters)

    df = pd.DataFrame(activities)

    if df.empty:
        logger.info(f"No bioactivity data found for {uniprot_id} ({target_chembl_id}) with readouts: {readouts}")
    else:
        logger.info(f"Retrieved {len(df)} bioactivities for {uniprot_id}")
    
    if not df.empty:
        df["uniprot_id"] = uniprot_id
        df["target_pref_name"] = target_query[0].get("pref_name")

    return df


@register_task("clean_bioactivities", 
               category='Bioactivity',
               description="Check and standardise bioactivities, and attach publication dates via CrossRef.")
def clean_bioactivities(config, data):
    import requests

    def fetch_crossref_date(doi):
        """Fetch full publication date (YYYY-MM-DD) from CrossRef given a DOI."""
        if not doi:
            return None
        url = f"https://api.crossref.org/works/{doi}"
        try:
            resp = requests.get(url, timeout=10)
            if resp.status_code != 200:
                return None
            msg = resp.json().get("message", {})
            date_parts = msg.get("issued", {}).get("date-parts", [])
            if date_parts and isinstance(date_parts[0], list):
                parts = date_parts[0]
                parts += [1] * (3 - len(parts))  
                return f"{parts[0]:04d}-{parts[1]:02d}-{parts[2]:02d}"
        except Exception as e:
            logger.warning(f"CrossRef lookup failed for DOI {doi}: {e}")
        return None

    uniprot_id = config.get("uniprot_id", "UNKNOWN")
    readout_priority = config.get("readout", ["IC50", "Ki", "EC50"])
    if isinstance(readout_priority, str):
        readout_priority = [readout_priority]

    if not isinstance(data, pd.DataFrame) or data.empty:
        logger.info(f"[{uniprot_id}] Received empty or invalid DataFrame.")
        return {"df": pd.DataFrame(), "readout": None}

    if "standard_type" not in data.columns:
        logger.info(f"❌ [{uniprot_id}] Missing 'standard_type' column.")
        return {"df": pd.DataFrame(), "readout": None}

    cleaned_frames = []
    for readout in tqdm(readout_priority, desc=f"[{uniprot_id}] Cleaning readouts", leave=False):
        selected_readout = readout.upper()
        df_readout = data[data["standard_type"].str.upper() == selected_readout]
        if df_readout.empty:
            logger.info(f"[{uniprot_id}] No data for readout '{selected_readout}'")
            continue

        if "standard_units" not in df_readout.columns:
            logger.info(f"❌ [{uniprot_id}] Missing 'standard_units' column.")
            continue

        df_readout = df_readout[df_readout["standard_units"].str.lower() == "nm"]

        if "standard_value" not in df_readout.columns or "molecule_chembl_id" not in df_readout.columns:
            logger.info(f"❌ [{uniprot_id}] Missing required columns.")
            continue

        df_readout = df_readout.dropna(subset=["standard_value", "molecule_chembl_id"])
        df_readout["standard_value"] = pd.to_numeric(df_readout["standard_value"], errors="coerce")
        df_readout = df_readout.dropna(subset=["standard_value"])

        if df_readout.empty:
            logger.info(f"[{uniprot_id}] No usable values for readout '{selected_readout}' after cleaning.")
            continue

        columns_to_keep = [
            "molecule_chembl_id", "target_chembl_id", "target_pref_name",
            "standard_type", "standard_relation", "standard_value", "standard_units",
            "assay_chembl_id", "assay_type", "assay_description", "document_chembl_id",
            "uniprot_id", "target_pref_name"
        ]
        df_readout = df_readout[[col for col in columns_to_keep if col in df_readout.columns]]
        df_readout["readout"] = selected_readout
        logger.info(f"[{uniprot_id}] Cleaned data for '{selected_readout}': shape = {df_readout.shape}")

        cleaned_frames.append(df_readout)

    if not cleaned_frames:
        logger.info(f"[{uniprot_id}] No readout data cleaned for any requested readouts.")
        return {"df": pd.DataFrame(), "readout": None}

    combined = pd.concat(cleaned_frames, ignore_index=True)

    # Because Chembl doesn't expose document dates, need to retrieve them.
    if "document_chembl_id" in combined.columns:
        document_ids = combined["document_chembl_id"].dropna().unique().tolist()
        if document_ids:
            logger.debug(f"[{uniprot_id}] Fetching document metadata for {len(document_ids)} entries...")
            document_client = new_client.document
            documents = document_client.filter(document_chembl_id__in=document_ids)
            doc_df = pd.DataFrame(documents)

            if not doc_df.empty and "doi" in doc_df.columns:
                logger.debug(f"[{uniprot_id}] Looking up publication dates via CrossRef...")
                doc_df["publication_date"] = doc_df["doi"].apply(fetch_crossref_date)

                merge_fields = ["document_chembl_id", "doi", "year", "publication_date"]
                doc_df = doc_df[[col for col in merge_fields if col in doc_df.columns]]
                doc_df = doc_df.drop_duplicates(subset="document_chembl_id")

                combined = pd.merge(combined, doc_df, on="document_chembl_id", how="left")
                logger.debug(f"[{uniprot_id}] Merged CrossRef publication dates into bioactivity data.")
            else:
                logger.warning(f"[{uniprot_id}] No DOIs available for CrossRef lookup.")

    return {"df": combined, "readout": None}



@register_task("retrieve_compound_smiles", 
               category='Bioactivity',
               description="Retrieve SMILES from downloaded compound data.")
def retrieve_compound_smiles(config, data):
    if not isinstance(data, dict) or "df" not in data:
        raise ValueError("Expected a dict with 'df' key containing a DataFrame.")

    input_df = data["df"]
    if input_df.empty:
        raise ValueError("Input DataFrame is empty.")

    logger.info(f"Fetching SMILES for {len(input_df)} entries...")

    molecule_ids = input_df["molecule_chembl_id"].dropna().unique().tolist()
    compounds_api = new_client.molecule

    def chunk_list(lst, chunk_size):
        for i in range(0, len(lst), chunk_size):
            yield lst[i:i + chunk_size]

    all_compounds = []
    for chunk in tqdm(chunk_list(molecule_ids, 25), total=math.ceil(len(molecule_ids) / 25), desc="Fetching compounds"):
        results = compounds_api.filter(molecule_chembl_id__in=chunk)
        all_compounds.extend(list(results))

    compound_df = pd.DataFrame.from_records(all_compounds)

    if compound_df.empty:
        raise ValueError("No compound data retrieved from ChEMBL.")

    # Extract canonical SMILES
    compound_df["smiles"] = compound_df["molecule_structures"].apply(
        lambda x: x.get("canonical_smiles") if isinstance(x, dict) else None
    )
    compound_df = compound_df.dropna(subset=["smiles"])
    compound_df = compound_df.drop_duplicates(subset="molecule_chembl_id")

    merged_df = pd.merge(input_df, compound_df[["molecule_chembl_id", "smiles"]], on="molecule_chembl_id", how="left")

    if merged_df["smiles"].isnull().all():
        raise ValueError("No SMILES could be attached to the input compounds.")

    logger.info(f"Attached SMILES to {merged_df['smiles'].notnull().sum()} entries.")

    return {"df": merged_df, "readout": data.get("readout")}

@register_task("annotate_bioactivity_pactivity", 
               category='Bioactivity',
               description="Compute p(readout)) and add to retrieval results.")
def annotate_bioactivity_pactivity(config, data):
    if not isinstance(data, dict) or "df" not in data:
        raise ValueError("Expected a dict with 'df' key containing a DataFrame.")

    df = data["df"]
    readout = data.get("readout")

    if df.empty:
        raise ValueError("Input DataFrame is empty.")

    if "standard_value" not in df.columns:
        raise ValueError("Missing 'standard_value' column in input data.")

    bioactivity_readouts = {
        "IC50", "EC50", "Ki", "Kd", "Potency",
        "AC50", "GI50", "MIC", "XC50", "Activity"
    }

    # Guess readout if not explicitly passed
    readout_guess = readout or df.get("readout", pd.Series()).mode().iloc[0] if "readout" in df.columns else None
    if readout_guess is None or readout_guess not in bioactivity_readouts:
        logger.info(f"[annotate_bioactivity_pactivity] Skipping pActivity. Not a bioactivity readout: {readout_guess}")
        return {"df": df, "readout": readout_guess}

    logger.info(f"Calculating pActivity for readout: {readout_guess}")

    def compute_pX(val):
        if pd.notnull(val) and val > 0:
            try:
                return 9 - math.log10(val)
            except:
                return None
        return None

    df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
    df["pActivity"] = df["standard_value"].apply(compute_pX)

    num_valid = df["pActivity"].notnull().sum()
    logger.info(f"Computed pActivity for {num_valid} entries.")

    return {"df": df, "readout": readout_guess}