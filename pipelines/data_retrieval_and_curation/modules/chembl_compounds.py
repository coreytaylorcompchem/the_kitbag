from pipeline.task_registry import register_task
from chembl_webresource_client.new_client import new_client

import pandas as pd
from tqdm import tqdm
import time

from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

@register_task("retrieve_compound_metadata_batch",
               category="Compound",
               description="Retrieve metadata for a batch of ChEMBL molecule IDs.")
def retrieve_compound_metadata_batch(config, data=None):
    molecule_ids = config.get("molecule_chembl_ids", [])
    if not molecule_ids:
        raise ValueError("No molecule_chembl_ids provided in config.")

    compounds_api = new_client.molecule

    def chunk_list(lst, chunk_size):
        for i in range(0, len(lst), chunk_size):
            yield lst[i:i + chunk_size]

    all_compounds = []
    for chunk in tqdm(chunk_list(molecule_ids, 50), desc="Fetching batch", leave=False):
        try:
            results = compounds_api.filter(molecule_chembl_id__in=chunk)
            all_compounds.extend(list(results))
        except Exception as e:
            logger.warning(f"Chunk failed: {e}")

    df = pd.DataFrame.from_records(all_compounds)
    if df.empty:
        logger.warning("No compounds retrieved for batch.")
        return pd.DataFrame()

    # Extract SMILES
    df["smiles"] = df["molecule_structures"].apply(
        lambda x: x.get("canonical_smiles") if isinstance(x, dict) else None
    )

    keep_cols = [
        "molecule_chembl_id", "pref_name", "molecule_type", "max_phase",
        "therapeutic_flag", "first_approval", "oral", "parenteral",
        "topical", "black_box_warning", "natural_product", "first_in_class",
        "smiles"
    ]
    return df[[col for col in keep_cols if col in df.columns]]
