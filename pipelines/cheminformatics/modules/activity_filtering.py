from pipeline.task_registry import register_task
from pipeline.logger import setup_logger
from pathlib import Path
import pandas as pd

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_task("activity_filtering", category="Filtering", description="Filter molecules based on activity cutoff")
def activity_filtering(config, data=None):
    input_file = config.get("input_file")
    if input_file is None:
        raise ValueError("No input_file specified in config")

    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file does not exist: {input_file}")

    df = pd.read_csv(input_file)
    if "smiles" not in df.columns or "standard_value" not in df.columns:
        raise ValueError("Input CSV must contain 'smiles' and 'standard_value' columns")

    activity_cutoff = config.get("activity", {}).get("value_cutoff_nM", 1000)

    # Drop NA and convert activity to numeric
    df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
    df = df.dropna(subset=["standard_value"])

    filtered_df = df[df["standard_value"] <= activity_cutoff].reset_index(drop=True)

    logger.debug(f"activity_filtering: {len(filtered_df)} / {len(df)} molecules retained (≤ {activity_cutoff} nM)")

    return (input_path.stem, filtered_df)
