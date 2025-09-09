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
    return pd.DataFrame.from_records(bioactivities)


@register_task("clean_bioactivities")
def clean_bioactivities(config, df):
    df = df.dropna()
    df = df[df["standard_units"] == "nM"]
    df = df.drop_duplicates("molecule_chembl_id")
    df.rename(columns={"standard_value": "IC50", "standard_units": "units"}, inplace=True)
    return df

@register_task("retrieve_compound_data")
def retrieve_compound_data(config, bio_df):
    compounds_api = new_client.molecule
    ids = list(bio_df["molecule_chembl_id"])
    # Request the necessary fields explicitly:
    compounds = list(
        compounds_api.filter(molecule_chembl_id__in=ids)
        .only("molecule_chembl_id", "molecule_structures")
    )
    df = pd.DataFrame.from_records(compounds)

    if "molecule_structures" not in df.columns:
        raise KeyError("The field 'molecule_structures' is missing from the retrieved compound data.")

    df["smiles"] = df["molecule_structures"].apply(lambda x: x.get("canonical_smiles") if x else None)
    df = df.dropna(subset=["smiles"])
    df = df.drop_duplicates("molecule_chembl_id")

    merged = pd.merge(bio_df, df[["molecule_chembl_id", "smiles"]], on="molecule_chembl_id")
    merged["pIC50"] = merged["IC50"].astype(float).apply(lambda x: 9 - math.log10(x))
    return merged
