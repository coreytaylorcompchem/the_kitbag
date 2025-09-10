from pipeline.task_registry import register_task
from chembl_webresource_client.new_client import new_client
import pandas as pd
import math
from tqdm import tqdm

@register_task("retrieve_chembl_bioactivities")
def retrieve_chembl_bioactivities(config, _):
    uniprot_id = config["uniprot_id"]
    readout = config["readout"]
    relation = config["relation"]
    assay_type = config["assay_type"]

    targets = new_client.target.get(target_components__accession=uniprot_id)
    targets = pd.DataFrame.from_records(targets)
    if targets.empty:
        raise ValueError("No targets found.")
    chembl_id = targets.iloc[0]["target_chembl_id"]

    bioactivities = new_client.activity.filter(
        target_chembl_id=chembl_id,
        type=readout,
        relation=relation,
        assay_type=assay_type,
    )

    bioactivities = list(tqdm(bioactivities))
    print(f"Retrieved {len(bioactivities)} bioactivities")
    return pd.DataFrame.from_records(bioactivities)


@register_task("clean_bioactivities")
def clean_bioactivities(config, df):
    # Print initial shape
    print(f"Cleaning bioactivities: input shape = {df.shape}")
    
    # Check required columns
    required_cols = ["standard_value", "standard_units", "molecule_chembl_id"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in bioactivity data: {missing_cols}")
    
    # Only keep rows with required fields not null
    df = df.dropna(subset=required_cols)
    
    # Focus on nanomolar values
    df = df[df["standard_units"] == "nM"]
    print(f"After filtering for 'nM': shape = {df.shape}")
    
    # Drop duplicates
    df = df.drop_duplicates("molecule_chembl_id")

    # Rename for downstream processing
    df.rename(columns={"standard_value": "IC50", "standard_units": "units"}, inplace=True)
    
    print(f"Final cleaned bioactivities: shape = {df.shape}")
    return df


@register_task("retrieve_compound_data")
def retrieve_compound_data(config, bio_df):

    compounds_api = new_client.molecule
    ids = list(bio_df["molecule_chembl_id"])

    # Handle batching to avoid overloading API
    def chunk_list(lst, n):
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    all_compounds = []
    for chunk in chunk_list(ids, 25):  # safer chunk size
        results = compounds_api.filter(molecule_chembl_id__in=chunk)
        records = list(results)
        if not records:
            print(f"No compound data returned for chunk: {chunk}")
        all_compounds.extend(records)

    df = pd.DataFrame.from_records(all_compounds)

    # Early exit if no data retrieved
    if df.empty:
        raise ValueError("No compound data retrieved from ChEMBL.")

    # Check for expected columns
    missing_cols = [col for col in ["molecule_chembl_id", "molecule_structures"] if col not in df.columns]
    if missing_cols:
        raise KeyError(f"Missing expected column(s) in compound data: {missing_cols}\nColumns present: {df.columns.tolist()}")

    # Extract SMILES
    try:
        df["smiles"] = df["molecule_structures"].apply(lambda x: x.get("canonical_smiles") if x else None)
    except Exception as e:
        print("Error extracting SMILES:", e)
    df = df.dropna(subset=["smiles"])
    df = df.drop_duplicates("molecule_chembl_id")

    if df.empty:
        raise ValueError("No compounds with valid SMILES found.")

    # Final merge
    merged = pd.merge(bio_df, df[["molecule_chembl_id", "smiles"]], on="molecule_chembl_id", how="inner")
    if merged.empty:
        raise ValueError("No overlap between bioactivity data and compound SMILES.")

    # drop any weird values
    merged["IC50"] = pd.to_numeric(merged["IC50"], errors="coerce")
    merged = merged.dropna(subset=["IC50"])
    merged["pIC50"] = merged["IC50"].apply(lambda x: 9 - math.log10(x) if x > 0 else None)

    output_df = merged[[
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
        "IC50", 
        "pIC50", 
        "units"
        ]]
    
    print(f"Final output shape: {output_df.shape}")

    return output_df
