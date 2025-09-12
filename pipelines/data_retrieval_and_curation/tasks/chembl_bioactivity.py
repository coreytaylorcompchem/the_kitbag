from pipeline.task_registry import register_task
from chembl_webresource_client.new_client import new_client
import pandas as pd
import math
from tqdm import tqdm

@register_task("retrieve_chembl_bioactivities", description="Retrieve bioactivity data from CHEMBL.")
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


@register_task("clean_bioactivities", description="Check and standardise bioactivities.")
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

        if "standard_units" not in df_readout.columns:
            print(f"❌ [{uniprot_id}] Missing 'standard_units' column.")
            continue

        df_readout = df_readout[df_readout["standard_units"].str.lower() == "nm"]

        # Drop NaN in compulsory columns
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


@register_task("retrieve_compound_smiles", description="Retrieve SMILES from downloaded compound data.")
def retrieve_compound_smiles(config, data):
    if not isinstance(data, dict) or "df" not in data:
        raise ValueError("Expected a dict with 'df' key containing a DataFrame.")

    input_df = data["df"]
    if input_df.empty:
        raise ValueError("Input DataFrame is empty.")

    print(f"[retrieve_compound_smiles] Fetching SMILES for {len(input_df)} entries...")

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

    print(f"[retrieve_compound_smiles] Attached SMILES to {merged_df['smiles'].notnull().sum()} entries.")

    return {"df": merged_df, "readout": data.get("readout")}

@register_task("annotate_bioactivity_pactivity", description="Compute p(readout)) and add to retrieval results.")
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
        print(f"[annotate_bioactivity_pactivity] Skipping pActivity. Not a bioactivity readout: {readout_guess}")
        return {"df": df, "readout": readout_guess}

    print(f"[annotate_bioactivity_pactivity] Calculating pActivity for readout: {readout_guess}")

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
    print(f"[annotate_bioactivity_pactivity] Computed pActivity for {num_valid} entries.")

    return {"df": df, "readout": readout_guess}
