import os
import math
import pandas as pd
from copy import deepcopy
from multiprocessing import Pool, cpu_count
from pathlib import Path
from tqdm import tqdm

from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

class ParallelWorkflowRunner:
    def __init__(
        self,
        workflow_func,
        config: dict,
        input_key: str = "uniprot_ids",
        output_key: str = "uniprot_id",
        output_dir: str = "outputs/parallel",
        filename_pattern: str = "{uniprot_id}_bioactivity.csv",
        combined_filename: str = "combined_output.csv",
        input_is_pair=False,
        use_multiprocessing: bool = True,
        reserve_cpus: int = 4
    ):
        self.workflow_func = workflow_func
        self.config = config
        self.input_key = input_key
        self.output_key = output_key
        self.output_dir = Path(config.get("output", {}).get("directory", output_dir))
        self.filename_pattern = filename_pattern
        self.combined_filename = config.get("output", {}).get("filename", combined_filename)
        self.use_multiprocessing = use_multiprocessing
        self.reserve_cpus = reserve_cpus
        self.input_is_pair = input_is_pair

        self.inputs = config.get(self.input_key, [])
        if not self.inputs:
            raise ValueError(f"No inputs provided in config under '{self.input_key}'.")

        self.output_dir.mkdir(parents=True, exist_ok=True)

    def _run_for_single(self, identifier):
        try:
            result = self.workflow_func(identifier)
        except Exception as e:
            logger.info(f"[{identifier}] Error: {e}")
            return pd.DataFrame()

        logger.debug(f"[{identifier}] Result type: {type(result)}")
        if isinstance(result, tuple):
            logger.debug(f"[{identifier}] Tuple length: {len(result)}")
            logger.debug(f"[{identifier}] Tuple content: {result}")
        elif isinstance(result, dict):
            logger.debug(f"[{identifier}] Dict keys: {list(result.keys())}")

        # Extract DataFrame and filename_key from known result structures
        if isinstance(result, tuple):
            # Defensive: ensure second element is DataFrame
            if len(result) >= 2 and isinstance(result[1], pd.DataFrame):
                filename_key, result_df = result[0], result[1]
            else:
                logger.info(f"[{identifier}] Tuple result does not contain expected DataFrame at index 1.")
                return pd.DataFrame()
        elif isinstance(result, pd.DataFrame):
            filename_key = identifier
            result_df = result
        elif isinstance(result, dict) and "df" in result and isinstance(result["df"], pd.DataFrame):
            filename_key = identifier
            result_df = result["df"]
        else:
            logger.info(f"[{identifier}] Unsupported result type: {type(result)}")
            return pd.DataFrame()

        if not isinstance(result_df, pd.DataFrame):
            logger.info(f"[{identifier}] result_df is not a DataFrame: {type(result_df)}")
            return pd.DataFrame()

        if result_df.empty:
            logger.debug(f"[{identifier}] No data returned (empty DataFrame).")
            return pd.DataFrame()

        # Sanitize filename_key string
        if not isinstance(filename_key, str):
            filename_key = str(filename_key)

        try:
            filename_path = Path(filename_key)
            if filename_path.suffix:
                filename_key = filename_path.stem
        except Exception as e:
            logger.warning(f"[{identifier}] Failed to parse filename_key '{filename_key}': {e}")
            filename_key = str(identifier)

        try:
            filename = self.filename_pattern.format(**{self.output_key: filename_key})
        except KeyError:
            filename = f"{filename_key}.csv"

        output_path = self.output_dir / filename
        try:
            result_df.to_csv(output_path, index=False)
            logger.debug(f"[{filename_key}] Saved: {output_path}")
        except Exception as e:
            logger.error(f"[{filename_key}] Failed to write CSV: {e}")
            return pd.DataFrame()

        logger.debug(f"[{identifier}] Returning DataFrame shape: {result_df.shape}")
        logger.debug(f"[{identifier}] DataFrame columns: {result_df.columns.tolist()}")
        return result_df

    def run(self):
        logger.info(f"Running workflow on {len(self.inputs)} inputs...")

        if self.use_multiprocessing and len(self.inputs) > 1:
            available_cpus = max(1, cpu_count() - self.reserve_cpus)
            num_workers = min(available_cpus, len(self.inputs))
            logger.info(f"Using {num_workers} worker(s) for parallel execution.")

            with Pool(processes=num_workers) as pool:
                results_iter = pool.imap_unordered(self._run_for_single, self.inputs)
                results = []
                for result in tqdm(results_iter, total=len(self.inputs), desc="Tasks"):
                    results.append(result)
        else:
            logger.info("Running in serial mode (no multiprocessing).")
            results = []
            for i in tqdm(self.inputs, desc="Tasks"):
                results.append(self._run_for_single(i))

        # Combine only non-empty DataFrames
        valid_dfs = [df for df in results if isinstance(df, pd.DataFrame) and not df.empty]
        if not valid_dfs:
            logger.info("No valid results to combine.")
            return pd.DataFrame()

        combined_df = pd.concat(valid_dfs, ignore_index=True)
        combined_path = self.output_dir / self.combined_filename
        try:
            combined_df.to_csv(combined_path, index=False)
            logger.info(f"Combined output saved to {combined_path}")
        except Exception as e:
            logger.error(f"Failed to write combined output CSV: {e}")
            return pd.DataFrame()

        return combined_df
