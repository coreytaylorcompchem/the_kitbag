from pipeline.task_registry import register_task
import pandas as pd
import re

@register_task("clean_adme_data")
def clean_adme_data(config, df=None):
    """
    Clean ADME data, agnostic of data source.
    Expects config dict with 'clean_adme_data' and 'cleaning' sections:
      - column_map: mapping of standard columns to data source columns
      - default_units: list of acceptable units for filtering
      - output_columns: list of columns to retain in final output
    """
    if df is None:
        raise ValueError("Input dataframe 'df' must be provided")

    # Unwrap nested df if needed
    if isinstance(df, dict) and "df" in df:
        df = df["df"]

    # Force conversion to DataFrame if dict
    if isinstance(df, dict):
        if all(not isinstance(v, (list, tuple, pd.Series)) for v in df.values()):
            df = pd.DataFrame([df])
        else:
            df = pd.DataFrame(df)

    # Load config sections
    cleaning_cfg = config.get("cleaning", {})
    task_cfg = config.get("clean_adme_data", {})
    default_units = cleaning_cfg.get("default_units", [])
    column_map = cleaning_cfg.get("column_map", {})
    output_columns = config.get("output_columns", [])

    print(f"[clean_adme_data] Loaded column_map: {column_map}")
    print(f"[clean_adme_data] Loaded default_units: {default_units}")
    print(f"[clean_adme_data] Output columns (if specified): {output_columns}")

    # Rename columns to standard names
    df = df.rename(columns=column_map)
    df = df.loc[:, ~df.columns.duplicated()]  # Remove duplicate columns if any
    print(f"[clean_adme_data] Columns after renaming: {df.columns.tolist()}")

    def normalize_unit(u):
        if isinstance(u, pd.Series):
            print(f"[normalize_unit] ERROR: Got Series instead of scalar:\n{u}")
            return ""
        if pd.isna(u):
            return ""
        u = str(u).lower().replace("µ", "u").replace("-", " ").replace("/", " ").replace(".", " ").strip()
        u = " ".join(u.split())  # normalize whitespace
        return u

    # Filter by allowed units
    if "standard_units" in df.columns:
        df["standard_units"] = df["standard_units"].astype(str)
        df["standard_units_norm"] = df["standard_units"].apply(normalize_unit)
        valid_units_norm = [normalize_unit(u) for u in default_units]

        print(f"[clean_adme_data] Normalized units in data: {df['standard_units_norm'].unique()}")
        print(f"[clean_adme_data] Normalized allowed units: {valid_units_norm}")

        def unit_match(unit_str):
            return any(allowed_unit in unit_str for allowed_unit in valid_units_norm)

        mask = df["standard_units_norm"].apply(unit_match)
        print(f"[clean_adme_data] Rows before filtering: {len(df)}, after filtering: {mask.sum()}")

        df = df[mask]

        # Assign final normalized column if needed in output
        if "normalized_units" in output_columns:
            df["normalized_units"] = df["standard_units_norm"]

        df = df.drop(columns=["standard_units_norm"])
    else:
        print("[clean_adme_data] WARNING: 'standard_units' column missing, skipping unit filtering.")

    # Retain only requested output columns
    if output_columns:
        missing_cols = [col for col in output_columns if col not in df.columns and col != "smiles"]
        if missing_cols:
            print(f"[clean_adme_data] WARNING: Missing output columns: {missing_cols}")
        df = df[[col for col in output_columns if col in df.columns]]

    readout = task_cfg.get("readout", None)

    return {
        "df": df,
        "readout": readout
    }
