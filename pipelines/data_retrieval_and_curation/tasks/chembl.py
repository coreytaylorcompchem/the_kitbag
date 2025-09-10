from pipeline.task_registry import register_task
from chembl_webresource_client.new_client import new_client
import pandas as pd
import math
from tqdm import tqdm

@register_task("retrieve_chembl_bioactivities")
def retrieve_chembl_bioactivities(config, data=None):
    uniprot_id = config.get("uniprot_id")
    assay_type = config.get("assay_type")
    relation = config.get("relation")
    readouts = config.get("readout", ["IC50"])
    
    if isinstance(readouts, str):
        readouts = [readouts]

    target = new_client.target
    activity = new_client.activity

    # Lookup target by uniprot
    target_query = target.filter(target_components__accession=uniprot_id)
    if not target_query:
        print(f"❌ No targets found for UniProt ID: {uniprot_id}")
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
    activities = activity.filter(**filters)

    df = pd.DataFrame(activities)

    if df.empty:
        print(f"No bioactivity data found for {uniprot_id} ({target_chembl_id}) with readouts: {readouts}")
    else:
        print(f"Retrieved {len(df)} bioactivities for {uniprot_id}")

    return df


@register_task("clean_bioactivities")
def clean_bioactivities(config, data):
    uniprot_id = config.get("uniprot_id", "UNKNOWN")
    readout_priority = config.get("readout", ["IC50", "Ki", "EC50"])
    if isinstance(readout_priority, str):
        readout_priority = [readout_priority]

    if not isinstance(data, pd.DataFrame) or data.empty:
        print(f"[{uniprot_id}] Received empty or invalid DataFrame.")
        return {"df": pd.DataFrame(), "readout": None}

    if "standard_type" not in data.columns:
        print(f"❌ [{uniprot_id}] Missing 'standard_type' column.")
        return {"df": pd.DataFrame(), "readout": None}

    cleaned_frames = []
    for readout in tqdm(readout_priority, desc=f"[{uniprot_id}] Cleaning readouts", leave=False):
        selected_readout = readout.upper()

        df_readout = data[data["standard_type"].str.upper() == selected_readout]
        if df_readout.empty:
            print(f"[{uniprot_id}] No data for readout '{selected_readout}'")
            continue

        # Filter for units nm
        if "standard_units" not in df_readout.columns:
            print(f"❌ [{uniprot_id}] Missing 'standard_units' column.")
            continue

        df_readout = df_readout[df_readout["standard_units"].str.lower() == "nm"]

        # Drop nulls in critical columns
        if "standard_value" not in df_readout.columns or "molecule_chembl_id" not in df_readout.columns:
            print(f"❌ [{uniprot_id}] Missing required columns.")
            continue

        df_readout = df_readout.dropna(subset=["standard_value", "molecule_chembl_id"])
        df_readout["standard_value"] = pd.to_numeric(df_readout["standard_value"], errors="coerce")
        df_readout = df_readout.dropna(subset=["standard_value"])

        if df_readout.empty:
            print(f"[{uniprot_id}] No usable values for readout '{selected_readout}' after cleaning.")
            continue

        columns_to_keep = [
            "molecule_chembl_id", "target_chembl_id", "target_pref_name",
            "standard_type", "standard_relation", "standard_value", "standard_units",
            "assay_chembl_id", "assay_type", "assay_description", "document_chembl_id"
        ]
        df_readout = df_readout[[col for col in columns_to_keep if col in df_readout.columns]]

        df_readout["readout"] = selected_readout

        print(f"[{uniprot_id}] Cleaned data for '{selected_readout}': shape = {df_readout.shape}")

        cleaned_frames.append(df_readout)

    if cleaned_frames:
        combined = pd.concat(cleaned_frames, ignore_index=True)
        return {"df": combined, "readout": None}
    else:
        print(f"[{uniprot_id}] No readout data cleaned for any requested readouts.")
        return {"df": pd.DataFrame(), "readout": None}


@register_task("retrieve_compound_data")
def retrieve_compound_data(config, data):
    if not isinstance(data, dict) or "df" not in data:
        raise ValueError(f"Expected dict with 'df' and 'readout', got {type(data)}")

    bio_df = data["df"]
    readout = data.get("readout")  # Can be None if multiple readouts

    if bio_df is None or bio_df.empty:
        raise ValueError("Bioactivity data is empty.")

    print(f"Received {len(bio_df)} entries to fetch compound data for.")
    
    compounds_api = new_client.molecule
    ids = list(bio_df["molecule_chembl_id"].unique())

    def chunk_list(lst, n):
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    all_compounds = []
    for chunk in tqdm(chunk_list(ids, 25), total=math.ceil(len(ids)/25), desc="Fetching compounds"):
        results = compounds_api.filter(molecule_chembl_id__in=chunk)
        records = list(results)
        all_compounds.extend(records)

    df = pd.DataFrame.from_records(all_compounds)

    if df.empty:
        raise ValueError("No compound data retrieved from ChEMBL.")

    df["smiles"] = df["molecule_structures"].apply(lambda x: x.get("canonical_smiles") if x else None)
    df = df.dropna(subset=["smiles"])
    df = df.drop_duplicates("molecule_chembl_id")

    merged = pd.merge(bio_df, df[["molecule_chembl_id", "smiles"]], on="molecule_chembl_id", how="inner")

    if merged.empty:
        raise ValueError("No overlap between bioactivity data and compound SMILES.")

    def compute_pX(row):
        val = row.get("standard_value")
        if pd.notnull(val) and val > 0:
            try:
                return 9 - math.log10(val)
            except:
                return None
        return None

    merged["standard_value"] = pd.to_numeric(merged["standard_value"], errors="coerce")
    merged["pActivity"] = merged.apply(compute_pX, axis=1)

    desired_cols = [
        "molecule_chembl_id", 
        "target_chembl_id", 
        "target_organism", 
        "target_pref_name",
        "molecule_pref_name",           
        "assay_description", 
        "bao_label", 
        "document_journal", 
        "document_year",
        "smiles", 
        "readout", 
        "standard_value", 
        "standard_units",
        "pActivity"
    ]

    output_df = merged[[col for col in desired_cols if col in merged.columns]]
    print(f"Final output shape: {output_df.shape}")

    return {"df": output_df}
