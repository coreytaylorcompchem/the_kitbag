import importlib
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors, Lipinski, Crippen
from rdkit.Chem import QED

from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED

from modules.utils.retrosynthetic_predictor_helpers import _process_building_block_chunk

from pipeline.task_registry import register_task

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_task(
    "load_building_blocks",
    category="Retrosenthetic predictor",
    description="Load commercial building blocks from a SMILES file."
)
def load_building_blocks(config, context):

    input_file = Path(config["input_file"])

    if not input_file.exists():
        raise FileNotFoundError(
            f"Building-block file not found: {input_file}"
        )

    smiles_col = config.get("smiles_col", "smiles")
    id_col = config.get("id_col", "mcule_id")

    records = []

    with open(input_file, "r") as f:
        for line_number, line in enumerate(f, start=1):

            line = line.strip()

            if not line:
                continue

            parts = line.split()

            if len(parts) < 2:
                logger.warning(
                    f"Skipping malformed line {line_number}: {line}"
                )
                continue

            records.append({
                smiles_col: parts[0],
                id_col: parts[1],
            })

    df = pd.DataFrame(records)

    if df.empty:
        raise ValueError(
            f"No building blocks loaded from {input_file}"
        )

    logger.info(
        f"Loaded {len(df):,} building blocks from {input_file}"
    )

    context["building_blocks"] = df

    return {
        "building_blocks": df
    }

@register_task(
    "process_building_blocks",
    category="SYNTHESIS",
    description="Canonicalise and calculate descriptors for commercial building blocks."
)
def process_building_blocks(config, context):

    df = context["building_blocks"].copy()

    smiles_col = config.get("smiles_col", "smiles")

    canonicalise = config.get("canonicalise_smiles", True)
    remove_invalid = config.get("remove_invalid_smiles", True)
    remove_duplicates = config.get("remove_duplicates", True)

    parallel_config = config.get("parallel", {})

    parallel_enabled = parallel_config.get("enabled", True)
    n_workers = parallel_config.get("n_workers", None)
    chunk_size = parallel_config.get("chunk_size", 10000)

    logger.info(
        f"Processing {len(df):,} building blocks"
    )

    # ------------------------------------------------------------------
    # Convert the input into chunks of ordinary Python dictionaries.
    # This avoids the very slow DataFrame.iterrows() loop.
    # ------------------------------------------------------------------

    records = df.to_dict(orient="records")

    chunks = (
        records[i:i + chunk_size]
        for i in range(0, len(records), chunk_size)
    )

    processed = []

    # ------------------------------------------------------------------
    # Parallel processing
    # ------------------------------------------------------------------

    if parallel_enabled:

        logger.info(
            f"Using multiprocessing with "
            f"{n_workers or 'default'} workers "
            f"and chunk size {chunk_size:,}"
        )

        with ProcessPoolExecutor(
            max_workers=n_workers
        ) as executor:

            pending = set()
            max_pending = (n_workers or 4) * 2

            with tqdm(
                total=len(df),
                desc="Processing building blocks"
            ) as progress:

                for chunk in chunks:

                    future = executor.submit(
                        _process_building_block_chunk,
                        chunk,
                        smiles_col,
                    )

                    pending.add(future)

                    # Keep the number of outstanding jobs bounded.
                    if len(pending) >= max_pending:

                        done, pending = wait(
                            pending,
                            return_when=FIRST_COMPLETED
                        )

                        for completed in done:

                            result = completed.result()

                            processed.extend(result)

                            progress.update(len(result))

                # Collect remaining jobs.
                while pending:

                    done, pending = wait(
                        pending,
                        return_when=FIRST_COMPLETED
                    )

                    for completed in done:

                        result = completed.result()

                        processed.extend(result)

                        progress.update(len(result))

    # ------------------------------------------------------------------
    # Serial fallback
    # ------------------------------------------------------------------

    else:

        logger.info(
            "Parallel processing disabled; using single process."
        )

        with tqdm(
            total=len(df),
            desc="Processing building blocks"
        ) as progress:

            for chunk in chunks:

                result = _process_building_block_chunk(
                    chunk,
                    smiles_col,
                )

                processed.extend(result)

                progress.update(len(result))

    # ------------------------------------------------------------------
    # Construct processed DataFrame
    # ------------------------------------------------------------------

    processed_df = pd.DataFrame(processed)

    # ------------------------------------------------------------------
    # Remove invalid molecules
    # ------------------------------------------------------------------

    if remove_invalid:

        n_invalid = (~processed_df["valid_smiles"]).sum()

        logger.info(
            f"Removing {n_invalid:,} invalid SMILES"
        )

        processed_df = processed_df[
            processed_df["valid_smiles"]
        ].copy()

    # ------------------------------------------------------------------
    # Remove duplicates
    # ------------------------------------------------------------------

    if remove_duplicates:

        before = len(processed_df)

        processed_df = processed_df.drop_duplicates(
            subset="inchikey"
        ).copy()

        logger.info(
            f"Removed {before - len(processed_df):,} duplicate molecules"
        )

    processed_df.reset_index(drop=True, inplace=True)

    logger.info(
        f"Processed building-block database: "
        f"{len(processed_df):,} unique molecules"
    )

    context["building_blocks"] = processed_df

    return {
        "building_blocks": processed_df
    }

@register_task(
    "save_building_block_database",
    category="SYNTHESIS",
    description="Save processed commercial building-block database."
)
def save_building_block_database(config, context):

    output_file = Path(config["output_file"])

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    df = context["building_blocks"]

    df.to_parquet(
        output_file,
        index=False
    )

    logger.info(
        f"Saved {len(df):,} building blocks to {output_file}"
    )

    return {
        "building_block_database": output_file
    }