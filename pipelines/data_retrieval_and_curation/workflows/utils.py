import pandas as pd

from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

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
        logger.info(f"[{readout}] ❌ Error: No DataFrame returned.")
        return pd.DataFrame()

    if not isinstance(df, pd.DataFrame):
        logger.info(f"[{readout}] ❌ Error: Returned object is not a DataFrame (got {type(df)}).")
        return pd.DataFrame()

    if df.empty:
        logger.info(f"[{readout}] ⚠️ Warning: DataFrame is empty after filtering. Possibly due to unit mismatch or filters.")
        return df  

    return df
