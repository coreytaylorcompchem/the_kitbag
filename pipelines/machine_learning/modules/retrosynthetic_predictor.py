import json
import math

from pathlib import Path

import pandas as pd
import numpy as np

from tqdm import tqdm

import matplotlib.pyplot as plt

from rdkit import Chem
# from rdkit.Chem import AllChem
# from rdkit.Chem import Descriptors, Lipinski, Crippen
# from rdkit.Chem import QED
from rdkit.Chem.Draw import rdMolDraw2D

from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED

from modules.utils.retrosynthetic_predictor_helpers import _process_building_block_chunk, _annotate_building_block_chunk, REACTION_HANDLE_SMARTS

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


@register_task(
    "validate_building_block_database",
    category="SYNTHESIS",
    description="Validate the processed commercial building-block database."
)
def validate_building_block_database(config, context):

    input_file = Path(
        config.get(
            "input_file",
            "outputs/commercial_bbs/mcule_building_blocks.parquet"
        )
    )

    output_file = Path(
        config.get(
            "output_file",
            "outputs/commercial_bbs/validation/validation_summary.json"
        )
    )

    if not input_file.exists():
        raise FileNotFoundError(
            f"Building-block database not found: {input_file}"
        )

    df = pd.read_parquet(input_file)

    logger.info(
        f"Validating building-block database: "
        f"{len(df):,} rows"
    )

    required_columns = [
        "smiles",
        "canonical_smiles",
        "inchikey",
        "valid_smiles",
        "molecular_weight",
        "logp",
        "tpsa",
        "hbd",
        "hba",
        "rotatable_bonds",
        "ring_count",
        "heavy_atom_count",
    ]

    missing_columns = [
        col for col in required_columns
        if col not in df.columns
    ]

    if missing_columns:
        raise ValueError(
            f"Missing required columns: {missing_columns}"
        )

    # ------------------------------------------------------------------
    # Basic integrity
    # ------------------------------------------------------------------

    n_rows = len(df)

    n_duplicate_inchikey = int(
        df["inchikey"].duplicated().sum()
    )

    n_duplicate_canonical_smiles = int(
        df["canonical_smiles"].duplicated().sum()
    )

    n_missing_smiles = int(
        df["smiles"].isna().sum()
    )

    n_missing_inchikey = int(
        df["inchikey"].isna().sum()
    )

    n_invalid_flag = int(
        (~df["valid_smiles"].astype(bool)).sum()
    )

    # ------------------------------------------------------------------
    # Multi-component molecules
    # ------------------------------------------------------------------

    is_multicomponent = (
        df["canonical_smiles"]
        .fillna("")
        .str.contains(r"\.", regex=True)
    )

    n_multicomponent = int(is_multicomponent.sum())

    # ------------------------------------------------------------------
    # Descriptor validation
    # ------------------------------------------------------------------

    descriptor_columns = [
        "molecular_weight",
        "logp",
        "tpsa",
        "hbd",
        "hba",
        "rotatable_bonds",
        "ring_count",
        "heavy_atom_count",
    ]

    descriptor_summary = {}

    for column in descriptor_columns:

        values = pd.to_numeric(
            df[column],
            errors="coerce"
        )

        finite = values[np.isfinite(values)]

        descriptor_summary[column] = {
            "missing": int(values.isna().sum()),
            "non_finite": int(
                (~np.isfinite(values.fillna(0))).sum()
            ),
            "min": float(finite.min()) if len(finite) else None,
            "max": float(finite.max()) if len(finite) else None,
            "mean": float(finite.mean()) if len(finite) else None,
            "median": float(finite.median()) if len(finite) else None,
            "q01": float(finite.quantile(0.01)) if len(finite) else None,
            "q05": float(finite.quantile(0.05)) if len(finite) else None,
            "q95": float(finite.quantile(0.95)) if len(finite) else None,
            "q99": float(finite.quantile(0.99)) if len(finite) else None,
        }

    # ------------------------------------------------------------------
    # Basic chemical sanity checks
    # ------------------------------------------------------------------

    very_small = df["molecular_weight"] < config.get(
        "very_small_mw",
        50.0
    )

    very_large = df["molecular_weight"] > config.get(
        "very_large_mw",
        800.0
    )

    zero_heavy_atoms = df["heavy_atom_count"] <= 0

    summary = {
        "input_file": str(input_file),
        "n_rows": n_rows,
        "n_unique_inchikey": int(df["inchikey"].nunique()),
        "n_duplicate_inchikey": n_duplicate_inchikey,
        "n_duplicate_canonical_smiles": n_duplicate_canonical_smiles,
        "n_missing_smiles": n_missing_smiles,
        "n_missing_inchikey": n_missing_inchikey,
        "n_invalid_smiles": n_invalid_flag,
        "n_multicomponent": n_multicomponent,
        "n_single_component": n_rows - n_multicomponent,
        "n_very_small": int(very_small.sum()),
        "n_very_large": int(very_large.sum()),
        "n_zero_heavy_atoms": int(zero_heavy_atoms.sum()),
        "descriptor_summary": descriptor_summary,
    }

    # ------------------------------------------------------------------
    # Save report
    # ------------------------------------------------------------------

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    with open(output_file, "w") as f:
        json.dump(
            summary,
            f,
            indent=2
        )

    logger.info(
        f"Validation complete. "
        f"{n_rows:,} rows, "
        f"{n_multicomponent:,} multi-component"
    )

    logger.info(
        f"Validation report saved to {output_file}"
    )

    return {
        "building_block_validation": summary,
        "building_block_validation_file": output_file,
    }

def plot_distribution_and_outliers(
    df,
    column,
    normal_min,
    normal_max,
    xlabel,
    title,
    output_dir,
    bins=100,
):
    """
    Generate two plots for a descriptor:

    1. Normal distribution within the specified range.
    2. Distribution of values outside that range.

    The ranges are for visualisation only and do not filter the
    underlying database.
    """

    values = pd.to_numeric(
        df[column],
        errors="coerce"
    ).dropna()

    normal = values[
        (values >= normal_min) &
        (values <= normal_max)
    ]

    outliers = values[
        (values < normal_min) |
        (values > normal_max)
    ]

    # --------------------------------------------------------------
    # Normal distribution
    # --------------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(10, 6)
    )

    if not normal.empty:
        ax.hist(
            normal,
            bins=bins,
        )

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Number of molecules")

    ax.set_title(
        f"{title} — normal range"
    )

    ax.set_xlim(
        normal_min,
        normal_max
    )

    ax.grid(
        alpha=0.2
    )

    fig.tight_layout()

    fig.savefig(
        output_dir / f"{column}_distribution.png",
        dpi=200
    )

    plt.close(fig)

    # --------------------------------------------------------------
    # Outlier distribution
    # --------------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(10, 6)
    )

    if not outliers.empty:

        # Don't let a handful of extreme values make this plot
        # unreadable either. Use the 1st–99th percentile for the
        # plotting range, but report the full outlier count.
        lower = outliers.quantile(0.01)
        upper = outliers.quantile(0.99)

        if lower == upper:
            lower = outliers.min()
            upper = outliers.max()

        if lower == upper:
            lower -= 1
            upper += 1

        ax.hist(
            outliers,
            bins=bins,
            range=(lower, upper),
        )

        ax.axvline(
            normal_min,
            linestyle="--",
            linewidth=1,
        )

        ax.axvline(
            normal_max,
            linestyle="--",
            linewidth=1,
        )

        ax.set_xlim(
            lower,
            upper
        )

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Number of molecules")

    ax.set_title(
        f"{title} — outliers "
        f"({len(outliers):,} molecules)"
    )

    ax.grid(
        alpha=0.2
    )

    fig.tight_layout()

    fig.savefig(
        output_dir / f"{column}_outliers.png",
        dpi=200
    )

    plt.close(fig)

    logger.info(
        f"{column}: "
        f"{len(normal):,} within normal range, "
        f"{len(outliers):,} outliers"
    )

@register_task(
    "plot_building_block_database",
    category="SYNTHESIS",
    description="Generate diagnostic plots and molecular depictions for the building-block database."
)
def plot_building_block_database(config, context):

    input_file = Path(
        config.get(
            "input_file",
            "outputs/commercial_bbs/mcule_building_blocks.parquet"
        )
    )

    output_dir = Path(
        config.get(
            "output_dir",
            "outputs/commercial_bbs/validation"
        )
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    if not input_file.exists():
        raise FileNotFoundError(
            f"Building-block database not found: {input_file}"
        )

    df = pd.read_parquet(input_file)

    logger.info(
        f"Generating plots for {len(df):,} building blocks"
    )

    # ------------------------------------------------------------------
    # Descriptor distributions
    # ------------------------------------------------------------------

    plot_distribution_and_outliers(
        df=df,
        column="molecular_weight",
        normal_min=config.get("mw_normal_min", 50),
        normal_max=config.get("mw_normal_max", 800),
        xlabel="Molecular weight",
        title="Building-block molecular weight",
        output_dir=output_dir,
    )

    plot_distribution_and_outliers(
        df=df,
        column="logp",
        normal_min=config.get("logp_normal_min", -5),
        normal_max=config.get("logp_normal_max", 10),
        xlabel="LogP",
        title="Building-block LogP",
        output_dir=output_dir,
    )

    plot_distribution_and_outliers(
        df=df,
        column="tpsa",
        normal_min=config.get("tpsa_normal_min", 0),
        normal_max=config.get("tpsa_normal_max", 250),
        xlabel="TPSA",
        title="Building-block TPSA",
        output_dir=output_dir,
    )

    plot_distribution_and_outliers(
        df=df,
        column="heavy_atom_count",
        normal_min=config.get("hac_normal_min", 3),
        normal_max=config.get("hac_normal_max", 60),
        xlabel="Heavy atom count",
        title="Building-block heavy atom count",
        output_dir=output_dir,
    )

    # ------------------------------------------------------------------
    # Multi-component vs single-component
    # ------------------------------------------------------------------

    is_multicomponent = (
        df["canonical_smiles"]
        .fillna("")
        .str.contains(r"\.")
    )

    counts = pd.Series(
        {
            "Single component": int((~is_multicomponent).sum()),
            "Multi-component": int(is_multicomponent.sum()),
        }
    )

    plt.figure(figsize=(8, 6))

    counts.plot(
        kind="bar"
    )

    plt.ylabel("Number of molecules")
    plt.title("Single-component vs multi-component structures")
    plt.xticks(rotation=0)
    plt.tight_layout()

    plt.savefig(
        output_dir / "multi_component_summary.png",
        dpi=200
    )

    plt.close()

    # ------------------------------------------------------------------
    # Molecular weight zoom
    #
    # This is useful because a handful of huge molecules can compress
    # the interesting part of the MW distribution.
    # ------------------------------------------------------------------

    mw_max = config.get(
        "mw_plot_max",
        1000
    )

    plt.figure(figsize=(10, 6))

    plt.hist(
        df.loc[
            df["molecular_weight"] <= mw_max,
            "molecular_weight"
        ].dropna(),
        bins=100
    )

    plt.xlabel("Molecular weight")
    plt.ylabel("Number of molecules")
    plt.title(
        f"Building-block molecular weight distribution "
        f"(MW ≤ {mw_max})"
    )

    plt.tight_layout()

    plt.savefig(
        output_dir / "molecular_weight_distribution_zoomed.png",
        dpi=200
    )

    plt.close()

    # ------------------------------------------------------------------
    # Representative / problematic molecules
    # ------------------------------------------------------------------

    problematic = []

    # Very small
    small = df[
        df["molecular_weight"] < config.get(
            "very_small_mw",
            50
        )
    ]

    problematic.append(
        ("Very small", small)
    )

    # Very large
    large = df[
        df["molecular_weight"] > config.get(
            "very_large_mw",
            800
        )
    ]

    problematic.append(
        ("Very large", large)
    )

    # Multi-component
    multi = df[is_multicomponent]

    problematic.append(
        ("Multi-component", multi)
    )

    # ------------------------------------------------------------------
    # Generate molecular depiction sheets
    # ------------------------------------------------------------------

    max_per_category = config.get(
        "max_molecules_per_category",
        25
    )

    for category, subset in problematic:

        if subset.empty:
            continue

        subset = subset.head(max_per_category)

        molecules = []
        legends = []

        for _, row in subset.iterrows():

            mol = Chem.MolFromSmiles(
                row["canonical_smiles"]
            )

            if mol is None:
                continue

            molecules.append(mol)

            legends.append(
                f"{row['inchikey']}\n"
                f"MW={row['molecular_weight']:.1f}"
            )

        if not molecules:
            continue

        safe_category = (
            category.lower()
            .replace(" ", "_")
            .replace("-", "_")
        )

        output_file = (
            output_dir /
            f"{safe_category}_molecules.svg"
        )

        # ------------------------------------------------------------------
        # RDKit SVG renderer
        #
        # This avoids the Cairo dependency required by MolsToGridImage().
        # ------------------------------------------------------------------

        mols_per_row = 5
        cell_width = 250
        cell_height = 250

        n_rows = math.ceil(
            len(molecules) / mols_per_row
        )

        drawer = rdMolDraw2D.MolDraw2DSVG(
            mols_per_row * cell_width,
            n_rows * cell_height,
            cell_width,
            cell_height,
        )

        drawer.DrawMolecules(
            molecules,
            legends=legends
        )

        drawer.FinishDrawing()

        svg = drawer.GetDrawingText()

        with open(output_file, "w") as f:
            f.write(svg)

        logger.info(
            f"Saved {len(molecules)} {category.lower()} "
            f"molecular depictions to {output_file}"
        )

    return {
        "building_block_plot_directory": output_dir,
    }

@register_task(
    "classify_building_blocks",
    category="SYNTHESIS",
    description="Classify commercial molecules and split the database into retrosynthesis and excluded sets."
)
def classify_building_blocks(config, context):

    input_file = Path(
        config.get(
            "input_file",
            "outputs/commercial_bbs/mcule_building_blocks.parquet"
        )
    )

    retrosynthesis_output = Path(
        config.get(
            "retrosynthesis_output",
            "outputs/commercial_bbs/mcule_retrosynthesis_building_blocks.parquet"
        )
    )

    excluded_output = Path(
        config.get(
            "excluded_output",
            "outputs/commercial_bbs/mcule_excluded_building_blocks.parquet"
        )
    )

    if not input_file.exists():
        raise FileNotFoundError(
            f"Building-block database not found: {input_file}"
        )

    df = pd.read_parquet(input_file)

    logger.info(
        f"Classifying {len(df):,} building blocks"
    )

    # ------------------------------------------------------------------
    # Configuration
    # ------------------------------------------------------------------

    very_small_mw = config.get(
        "very_small_mw",
        50.0
    )

    very_large_mw = config.get(
        "very_large_mw",
        800.0
    )

    min_heavy_atoms = config.get(
        "min_heavy_atoms",
        3
    )

    max_heavy_atoms = config.get(
        "max_heavy_atoms",
        60
    )

    allowed_elements = set(
        config.get(
            "allowed_elements",
            [
                "H",
                "B",
                "C",
                "N",
                "O",
                "F",
                "P",
                "S",
                "Cl",
                "Br",
                "I",
                "Si",
            ]
        )
    )

    # ------------------------------------------------------------------
    # Start with default classification
    # ------------------------------------------------------------------

    df["building_block_category"] = (
        "retrosynthesis_building_block"
    )

    df["retrosynthesis_usable"] = True

    # ------------------------------------------------------------------
    # Multi-component structures
    # ------------------------------------------------------------------

    multi_component = (
        df["canonical_smiles"]
        .fillna("")
        .str.contains(r"\.")
    )

    df.loc[
        multi_component,
        "building_block_category"
    ] = "multi_component"

    df.loc[
        multi_component,
        "retrosynthesis_usable"
    ] = False

    # ------------------------------------------------------------------
    # Very small molecules
    # ------------------------------------------------------------------

    very_small = (
        df["molecular_weight"] < very_small_mw
    )

    df.loc[
        very_small & df["retrosynthesis_usable"],
        "building_block_category"
    ] = "very_small"

    df.loc[
        very_small,
        "retrosynthesis_usable"
    ] = False

    # ------------------------------------------------------------------
    # Very large molecules
    # ------------------------------------------------------------------

    very_large = (
        df["molecular_weight"] > very_large_mw
    )

    df.loc[
        very_large & df["retrosynthesis_usable"],
        "building_block_category"
    ] = "very_large"

    df.loc[
        very_large,
        "retrosynthesis_usable"
    ] = False

    # ------------------------------------------------------------------
    # Heavy atom count
    # ------------------------------------------------------------------

    too_few_atoms = (
        df["heavy_atom_count"] < min_heavy_atoms
    )

    df.loc[
        too_few_atoms & df["retrosynthesis_usable"],
        "building_block_category"
    ] = "too_few_heavy_atoms"

    df.loc[
        too_few_atoms,
        "retrosynthesis_usable"
    ] = False

    too_many_atoms = (
        df["heavy_atom_count"] > max_heavy_atoms
    )

    df.loc[
        too_many_atoms & df["retrosynthesis_usable"],
        "building_block_category"
    ] = "too_many_heavy_atoms"

    df.loc[
        too_many_atoms,
        "retrosynthesis_usable"
    ] = False

    # ------------------------------------------------------------------
    # Unusual elements
    # ------------------------------------------------------------------

    unusual_elements = ~df["elements"].apply(
        lambda elements: set(elements).issubset(allowed_elements)
    )

    df.loc[
        unusual_elements & df["retrosynthesis_usable"],
        "building_block_category"
    ] = "unusual_elements"

    df.loc[
        unusual_elements,
        "retrosynthesis_usable"
    ] = False

    # ------------------------------------------------------------------
    # Save the two datasets
    # ------------------------------------------------------------------

    retrosynthesis_df = df[
        df["retrosynthesis_usable"]
    ].copy()

    excluded_df = df[
        ~df["retrosynthesis_usable"]
    ].copy()

    retrosynthesis_output.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    excluded_output.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    retrosynthesis_df.to_parquet(
        retrosynthesis_output,
        index=False
    )

    excluded_df.to_parquet(
        excluded_output,
        index=False
    )

    # ------------------------------------------------------------------
    # Report
    # ------------------------------------------------------------------

    category_counts = (
        df["building_block_category"]
        .value_counts()
        .to_dict()
    )

    logger.info(
        f"Retrosynthesis-usable: "
        f"{len(retrosynthesis_df):,}"
    )

    logger.info(
        f"Excluded: "
        f"{len(excluded_df):,}"
    )

    logger.info(
        "Building-block categories:"
    )

    for category, count in category_counts.items():
        logger.info(
            f"  {category}: {count:,}"
        )

    logger.info(
        f"Saved retrosynthesis database to "
        f"{retrosynthesis_output}"
    )

    logger.info(
        f"Saved excluded database to "
        f"{excluded_output}"
    )

    return {
        "building_blocks": retrosynthesis_df,
        "retrosynthesis_building_blocks": retrosynthesis_df,
        "excluded_building_blocks": excluded_df,
        "building_block_category_counts": category_counts,
        "retrosynthesis_building_block_database": retrosynthesis_output,
        "excluded_building_block_database": excluded_output,
    }

@register_task(
    "annotate_building_block_reactivity",
    category="SYNTHESIS",
    description="Annotate commercial building blocks with retrosynthetically useful reaction handles."
)
def annotate_building_block_reactivity(config, context):

    input_file = Path(
        config.get(
            "input_file",
            "outputs/commercial_bbs/mcule_retrosynthesis_building_blocks.parquet"
        )
    )

    output_file = Path(
        config.get(
            "output_file",
            "outputs/commercial_bbs/mcule_retrosynthesis_building_blocks_reactivity.parquet"
        )
    )

    if not input_file.exists():
        raise FileNotFoundError(
            f"Retrosynthesis building-block database not found: {input_file}"
        )

    df = pd.read_parquet(input_file)

    logger.info(
        f"Annotating reaction handles for {len(df):,} building blocks"
    )

    required_columns = [
        "canonical_smiles",
        "inchikey",
    ]

    missing_columns = [
        col for col in required_columns
        if col not in df.columns
    ]

    if missing_columns:
        raise ValueError(
            f"Missing required columns: {missing_columns}"
        )

    parallel_config = config.get("parallel", {})

    parallel_enabled = parallel_config.get(
        "enabled",
        True
    )

    n_workers = parallel_config.get(
        "n_workers",
        None
    )

    chunk_size = parallel_config.get(
        "chunk_size",
        10000
    )

    records = df.to_dict(orient="records")

    chunks = (
        records[i:i + chunk_size]
        for i in range(
            0,
            len(records),
            chunk_size
        )
    )

    annotated = []

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

            max_pending = (
                (n_workers or 4) * 2
            )

            with tqdm(
                total=len(df),
                desc="Annotating reaction handles"
            ) as progress:

                for chunk in chunks:

                    future = executor.submit(
                        _annotate_building_block_chunk,
                        chunk,
                    )

                    pending.add(future)

                    if len(pending) >= max_pending:

                        done, pending = wait(
                            pending,
                            return_when=FIRST_COMPLETED,
                        )

                        for completed in done:

                            result = completed.result()

                            annotated.extend(result)

                            progress.update(
                                len(result)
                            )

                while pending:

                    done, pending = wait(
                        pending,
                        return_when=FIRST_COMPLETED,
                    )

                    for completed in done:

                        result = completed.result()

                        annotated.extend(result)

                        progress.update(
                            len(result)
                        )

    else:

        logger.info(
            "Parallel processing disabled; "
            "using single process."
        )

        with tqdm(
            total=len(df),
            desc="Annotating reaction handles"
        ) as progress:

            for chunk in chunks:

                result = _annotate_building_block_chunk(
                    chunk
                )

                annotated.extend(result)

                progress.update(
                    len(result)
                )

    annotated_df = pd.DataFrame(
        annotated
    )

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    annotated_df.to_parquet(
        output_file,
        index=False
    )

    handle_counts = {}

    for handle in REACTION_HANDLE_SMARTS:

        count = int(
            annotated_df["reaction_handles"]
            .fillna("")
            .str.contains(
                rf"(^|;){handle}(;|$)",
                regex=True
            )
            .sum()
        )

        handle_counts[handle] = count

    logger.info(
        f"Building blocks with at least one "
        f"reaction handle: "
        f"{int(annotated_df['has_reaction_handle'].sum()):,}"
    )

    logger.info(
        f"Building blocks with no recognised "
        f"reaction handle: "
        f"{int((~annotated_df['has_reaction_handle']).sum()):,}"
    )

    logger.info(
        "Reaction-handle counts:"
    )

    for handle, count in handle_counts.items():

        logger.info(
            f"  {handle}: {count:,}"
        )

    logger.info(
        f"Saved annotated building-block database "
        f"to {output_file}"
    )

    context[
        "retrosynthesis_building_blocks"
    ] = annotated_df

    return {
        "retrosynthesis_building_blocks": annotated_df,
        "annotated_building_block_database": output_file,
        "reaction_handle_counts": handle_counts,
    }

@register_task(
    "build_building_block_search_index",
    category="SYNTHESIS",
    description="Build searchable reaction-handle indexes for commercial building blocks."
)
def build_building_block_search_index(config, context):

    input_file = Path(
        config.get(
            "input_file",
            "outputs/commercial_bbs/mcule_retrosynthesis_building_blocks_reactivity.parquet"
        )
    )

    output_dir = Path(
        config.get(
            "output_dir",
            "outputs/commercial_bbs/search_index"
        )
    )

    if not input_file.exists():
        raise FileNotFoundError(
            f"Annotated building-block database not found: {input_file}"
        )

    df = pd.read_parquet(input_file)

    logger.info(
        f"Building search index from {len(df):,} building blocks"
    )

    required_columns = [
        "canonical_smiles",
        "inchikey",
        "reaction_handles",
    ]

    missing_columns = [
        col for col in required_columns
        if col not in df.columns
    ]

    if missing_columns:
        raise ValueError(
            f"Missing required columns: {missing_columns}"
        )

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    # ------------------------------------------------------------------
    # Columns retained in the search index
    # ------------------------------------------------------------------

    index_columns = [
        "inchikey",
        "canonical_smiles",
    ]

    optional_columns = [
        "smiles",
        "mcule_id",
        "molecular_weight",
        "logp",
        "tpsa",
        "hbd",
        "hba",
        "rotatable_bonds",
        "ring_count",
        "heavy_atom_count",
    ]

    index_columns.extend(
        [
            column
            for column in optional_columns
            if column in df.columns
        ]
    )

    # ------------------------------------------------------------------
    # Expand reaction handles
    # ------------------------------------------------------------------

    index_df = df[
        index_columns + ["reaction_handles"]
    ].copy()

    index_df["reaction_handle"] = (
        index_df["reaction_handles"]
        .fillna("")
        .str.split(";")
    )

    index_df = index_df.explode(
        "reaction_handle"
    )

    index_df["reaction_handle"] = (
        index_df["reaction_handle"]
        .astype(str)
        .str.strip()
    )

    # Remove rows where no reaction handle exists.

    index_df = index_df[
        index_df["reaction_handle"] != ""
    ].copy()

    # The original combined string is no longer needed.

    index_df.drop(
        columns=["reaction_handles"],
        inplace=True
    )

    # Avoid accidental duplicate handle/molecule pairs.

    index_df.drop_duplicates(
        subset=[
            "reaction_handle",
            "inchikey",
        ],
        inplace=True,
    )

    index_df.reset_index(
        drop=True,
        inplace=True
    )

    # ------------------------------------------------------------------
    # Save master handle index
    # ------------------------------------------------------------------

    master_index_file = (
        output_dir /
        "building_blocks_by_reaction_handle.parquet"
    )

    index_df.to_parquet(
        master_index_file,
        index=False
    )

    logger.info(
        f"Saved master reaction-handle index "
        f"({len(index_df):,} entries) "
        f"to {master_index_file}"
    )

    # ------------------------------------------------------------------
    # Save one index per reaction handle
    # ------------------------------------------------------------------

    handle_files = {}

    for handle, handle_df in index_df.groupby(
        "reaction_handle",
        sort=True
    ):

        handle_df = handle_df.drop(
            columns=["reaction_handle"]
        ).reset_index(drop=True)

        handle_file = (
            output_dir /
            f"{handle}.parquet"
        )

        handle_df.to_parquet(
            handle_file,
            index=False
        )

        handle_files[handle] = str(
            handle_file
        )

        logger.info(
            f"  {handle}: "
            f"{len(handle_df):,} building blocks"
        )

    # ------------------------------------------------------------------
    # Save metadata
    # ------------------------------------------------------------------

    metadata = {
        "input_file": str(input_file),
        "n_building_blocks": int(len(df)),
        "n_index_entries": int(len(index_df)),
        "n_reaction_handles": int(
            index_df["reaction_handle"].nunique()
        ),
        "reaction_handles": sorted(
            index_df["reaction_handle"].unique().tolist()
        ),
        "handle_files": handle_files,
    }

    metadata_file = (
        output_dir /
        "search_index_metadata.json"
    )

    import json

    with open(
        metadata_file,
        "w"
    ) as f:
        json.dump(
            metadata,
            f,
            indent=2
        )

    logger.info(
        f"Saved search-index metadata "
        f"to {metadata_file}"
    )

    return {
        "building_block_search_index": index_df,
        "building_block_search_index_file": master_index_file,
        "building_block_search_index_dir": output_dir,
        "building_block_search_index_metadata": metadata_file,
    }