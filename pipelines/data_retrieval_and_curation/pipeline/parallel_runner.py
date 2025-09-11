import os
import math
import pandas as pd
from copy import deepcopy
from multiprocessing import Pool, cpu_count
from pathlib import Path
from tqdm import tqdm

class ParallelWorkflowRunner:
    def __init__(
        self,
        workflow_func,
        config: dict,
        input_key: str = "uniprot_ids",
        output_key: str = "uniprot_id",  # key used per instance
        output_dir: str = "outputs/parallel",
        filename_pattern: str = "{uniprot_id}_bioactivity.csv",
        combined_filename: str = "combined_output.csv",
        input_is_pair=False,  # <-- new flag
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
        if self.input_is_pair:
            # identifier is assumed to be a tuple: (input_value, config)
            try:
                result = self.workflow_func(identifier)
            except Exception as e:
                print(f"❌ [{identifier}] Error: {e}")
                return pd.DataFrame()
        else:
            local_config = deepcopy(self.config)
            local_config[self.output_key] = identifier
            try:
                result = self.workflow_func(local_config)
            except Exception as e:
                print(f"❌ [{identifier}] Error: {e}")
                return pd.DataFrame()

        if isinstance(result, dict) and "df" in result:
            result_df = result["df"]
        elif isinstance(result, pd.DataFrame):
            result_df = result
        else:
            print(f"[{identifier}] Unsupported result type: {type(result)}")
            return pd.DataFrame()

        if result_df.empty:
            print(f"[{identifier}] No data returned.")
            return pd.DataFrame()

        # filename pattern formatting for the input, works for both str or tuple keys
        if isinstance(identifier, tuple):
            # Use first element of tuple for filename formatting (e.g. the readout string)
            filename_key = identifier[0]
        else:
            filename_key = identifier

        filename = self.filename_pattern.format(**{self.output_key: filename_key})
        output_path = self.output_dir / filename
        result_df.to_csv(output_path, index=False)
        print(f"[{filename_key}] Saved: {output_path}")
        return result_df

    def run(self):
        print(f"Running workflow on {len(self.inputs)} inputs...")

        if self.use_multiprocessing and len(self.inputs) > 1:
            available_cpus = max(1, cpu_count() - self.reserve_cpus)
            num_workers = min(available_cpus, len(self.inputs))
            print(f"Using {num_workers} worker(s) for parallel execution.")

            with Pool(processes=num_workers) as pool:
                results_iter = pool.imap_unordered(self._run_for_single, self.inputs)
                results = []
                for result in tqdm(results_iter, total=len(self.inputs), desc="Tasks"):
                    results.append(result)

        else:
            print("Running in serial mode (no multiprocessing).")
            results = []
            for i in tqdm(self.inputs, desc="Tasks"):
                results.append(self._run_for_single(i))

        valid_dfs = [df for df in results if isinstance(df, pd.DataFrame) and not df.empty]
        if not valid_dfs:
            print("❌ No valid results to combine.")
            return pd.DataFrame()

        combined_df = pd.concat(valid_dfs, ignore_index=True)
        combined_path = self.output_dir / self.combined_filename
        combined_df.to_csv(combined_path, index=False)
        print(f"\n✅ Combined output saved to {combined_path}")

        return combined_df