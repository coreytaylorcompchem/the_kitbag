import pandas as pd

def process_readout_dataframe(readout: str, df: pd.DataFrame) -> pd.DataFrame:
    """
    Gracefully handle empty or invalid DataFrames for a given readout.
    
    Args:
        readout (str): Name of the readout.
        df (pd.DataFrame): The result DataFrame for the readout.
        
    Returns:
        pd.DataFrame: The cleaned or empty DataFrame.
    """
    if df is None:
        print(f"[{readout}] ❌ Error: No DataFrame returned.")
        return pd.DataFrame()

    if not isinstance(df, pd.DataFrame):
        print(f"[{readout}] ❌ Error: Returned object is not a DataFrame (got {type(df)}).")
        return pd.DataFrame()

    if df.empty:
        print(f"[{readout}] ⚠️ Warning: DataFrame is empty after filtering. Possibly due to unit mismatch or filters.")
        return df  # Returning empty DataFrame is fine — it's expected sometimes.

    return df
