import gc
import pandas as pd
import csv
import shutil
import inspect

from tqdm import tqdm
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
import pyarrow.parquet as pq
from sklearn.decomposition import IncrementalPCA

from pipeline.task_registry import get_task, list_tasks
from pipeline.parallel_runner import ParallelWorkflowRunner
from pipeline.utils.config_helpers import update_config_input_file
from pipeline.logger import setup_logger
from workflows import register_workflow

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def stream_parquet_rowgroups(parquet_path):
    """Yield (row_group_index, pandas.DataFrame) from a Parquet file."""
    pf = pq.ParquetFile(parquet_path)
    for i in range(pf.num_row_groups):
        df_chunk = pf.read_row_group(i).to_pandas()
        yield i, df_chunk

def load_input_file(input_path, max_rows=None):
    """
    Loads CSV, TSV, or Parquet depending on file extension.
    Optionally limits to `max_rows` for testing or memory control.
    """
    input_path = Path(input_path)
    suffix = input_path.suffix.lower()

    logger.info(f"[dynamic_task_runner] Loading input file: {input_path}")

    if suffix == ".csv":
        df = pd.read_csv(input_path, nrows=max_rows)
    elif suffix == ".tsv":
        df = pd.read_csv(input_path, sep="\t", nrows=max_rows)
    elif suffix == ".parquet":
        try:
            pf = pq.ParquetFile(input_path)

            if max_rows:
                df = pf.read_row_groups([0], columns=None).to_pandas()
                df = df.head(max_rows)
            else:
                row_groups = []
                for rg in pf.iter_batches():
                    row_groups.append(rg.to_pandas())
                df = pd.concat(row_groups, ignore_index=True)
        except ImportError:
            df = pd.read_parquet(input_path)
    else:
        raise ValueError(f"Unsupported input format: {suffix}")

    logger.info(f"[dynamic_task_runner] Loaded dataframe shape: {df.shape}")
    return df


class ChunkWrapper:
    def __init__(self, task_list, config):
        self.task_list = task_list
        self.config = config

    def __call__(self, chunk_file):
        return run_task_chain_on_chunk(chunk_file, self.task_list, self.config)

def clean_smiles(smiles_series):
    # remove tabs/newlines, strip spaces, remove stray quotes, and clean up just enough
    s = (
        smiles_series.astype(str)
        .str.replace(r"[\t\r\n]", "", regex=True)
        .str.strip()
    )
    # Remove spaces inside slashes but leave other spaces
    s = s.str.replace(r"\s*/\s*", "/", regex=True)
    return s

def run_task_chain_on_chunk(chunk_file: str, task_list, config):
    df = pd.read_csv(chunk_file, delimiter=",", quotechar='"', engine="python")
    current_file = chunk_file
    current_data = {"df": df}
    mols_out = []

    for task_name in task_list:
        task = get_task(task_name)
        if task is None:
            raise ValueError(f"Task '{task_name}' not found")

        task_config = dict(config)
        task_config["input_file"] = current_file

        logger.debug(f"Running {task_name} on chunk {chunk_file}")
        result = task(task_config, data=current_data)

        # Normalise result to a dict
        if isinstance(result, tuple):
            logger.debug(f"Task '{task_name}' returned tuple of length {len(result)}")

            if len(result) == 3:
                _, df_result, mols_result = result
                result = {"df": df_result, "mols": mols_result}
            elif len(result) == 2:
                _, df_result = result
                result = {"df": df_result, "mols": []}
            else:
                logger.debug(f"Task '{task_name}' returned unexpected tuple length {len(result)}")
                return (Path(chunk_file).stem, pd.DataFrame(), [])

        elif isinstance(result, dict):
            if "df" not in result or result["df"] is None or result["df"].empty:
                logger.debug(f"Task '{task_name}' returned empty or invalid dataframe")
                return (Path(chunk_file).stem, pd.DataFrame(), [])
        else:
            logger.error(f"Task '{task_name}' returned unsupported result type: {type(result)}")
            return (Path(chunk_file).stem, pd.DataFrame(), [])

        # At this point, result is guaranteed to be a valid dict with a DataFrame
        df = result["df"]
        mols = result.get("mols", [])
        if mols:
            mols_out.extend(mols)

        # Save intermediate result (with proper quoting to protect SMILES)
        temp_output = Path(chunk_file).with_name(f"{Path(chunk_file).stem}_{task_name}.csv")
        df.to_csv(temp_output, index=False, quoting=csv.QUOTE_ALL)
        current_file = str(temp_output)
        current_data = result

        gc.collect()

    return (Path(chunk_file).stem, df, mols_out)

@register_workflow("dynamic_task_runner", description="Run a task list in parallel on input chunks and output CSV/SDF/SMI")
def dynamic_task_runner(config):
    import inspect

    task_list = config.get("workflow", [])
    if not task_list:
        raise ValueError("No 'workflow' list of tasks specified in config")

    input_file = config.get("input_file")
    if not input_file:
        def task_requires_input_file(task_name):
            task_func = get_task(task_name)
            if not task_func:
                raise ValueError(f"Task '{task_name}' not found")
            try:
                sig = inspect.signature(task_func)
                param = sig.parameters.get("data")
                return param is not None and param.default == inspect.Parameter.empty
            except Exception as e:
                logger.warning(f"Could not inspect task '{task_name}': {e}")
                return True  # Default to True if inspection fails

        if any(task_requires_input_file(task_name) for task_name in task_list):
            raise ValueError("Missing 'input_file' in config, but required by one or more tasks.")
        else:
            logger.info("No 'input_file' needed for this workflow.")
            current_data = {}
            for task_name in task_list:
                task_func = get_task(task_name)
                logger.info(f"Running task '{task_name}' (no input_file needed)...")
                result = task_func(config, data=current_data)
                if isinstance(result, dict):
                    current_data = result
                else:
                    current_data = {"df": result}
            return current_data  # Done here for no-input workflows

    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    chunk_size = config.get("chunk_size", None)

    if not chunk_size or chunk_size <= 0:
        # === Run all tasks sequentially on the full dataset ===

        logger.info("No chunk_size specified or <= 0: Running all tasks sequentially on full dataset")

        df = load_input_file(input_file)
        if "smiles" not in df.columns:
            raise ValueError("Input file must contain 'smiles' column")
        df["smiles"] = clean_smiles(df["smiles"])
        df = df[df["smiles"].str.len() > 0]

        current_data = {"df": df}

        for task_name in task_list:
            task_func = get_task(task_name)
            logger.info(f"Running task '{task_name}' sequentially on full dataset...")
            result = task_func(config, data=current_data)

            if isinstance(result, tuple):
                if len(result) == 3:
                    _, df_result, mols_result = result
                    current_data = {"df": df_result, "mols": mols_result}
                elif len(result) == 2:
                    _, df_result = result
                    current_data = {"df": df_result, "mols": []}
                else:
                    raise ValueError(f"Unexpected tuple result length from task '{task_name}'")
            elif isinstance(result, dict):
                current_data = result
            else:
                current_data = {"df": result}

        combined_df = current_data.get("df", pd.DataFrame())
        all_mols = current_data.get("mols", [])

        output_dir = Path(config.get("output", {}).get("directory", "outputs/flexible_parallel"))
        output_dir.mkdir(parents=True, exist_ok=True)

        final_csv = output_dir / config.get("output", {}).get("filename", "final_output.csv")

        try:
            combined_df.to_csv(
                final_csv,
                index=False,
                quoting=csv.QUOTE_ALL,
                quotechar='"',
                doublequote=True,
                escapechar=None,
                lineterminator='\n'
            )
            logger.info(f"Final CSV written: {final_csv}")
        except Exception as e:
            logger.error(f"Failed to write final CSV: {e}")

        if "smiles" in combined_df.columns:
            smi_file = final_csv.with_suffix(".smi")
            try:
                with open(smi_file, "w") as f:
                    for smi in combined_df["smiles"].dropna():
                        f.write(smi + "\n")
                logger.info(f"SMILES file written: {smi_file}")
            except Exception as e:
                logger.error(f"Failed to write SMILES file: {e}")

        if all_mols:
            sdf_file = final_csv.with_suffix(".sdf")
            try:
                writer = Chem.SDWriter(str(sdf_file))
                for mol in all_mols:
                    writer.write(mol)
                writer.close()
                logger.info(f"SDF file written: {sdf_file}")
            except Exception as e:
                logger.error(f"Failed to write SDF file: {e}")
        else:
            logger.debug("No RDKit mols to write SDF.")

        return {"df": combined_df}

    # === ELSE do chunking + parallelism ===

    df = load_input_file(input_file)
    if "smiles" not in df.columns:
        raise ValueError("Input file must contain 'smiles' column")

    df["smiles"] = clean_smiles(df["smiles"])
    df = df[df["smiles"].str.len() > 0]

    smiles_list = df["smiles"].dropna().unique().tolist()

    output_dir = Path(config.get("output", {}).get("directory", "outputs/flexible_parallel"))
    output_dir.mkdir(parents=True, exist_ok=True)

    chunk_dir = output_dir / "chunks"
    chunk_dir.mkdir(parents=True, exist_ok=True)

    chunk_files = []
    for i in range(0, len(df), chunk_size):
        chunk_df = df.iloc[i: i + chunk_size].copy()
        chunk_df["smiles"] = clean_smiles(chunk_df["smiles"])
        chunk_file = chunk_dir / f"chunk_{i // chunk_size}.csv"
        chunk_df.to_csv(
            chunk_file,
            index=False,
            quoting=csv.QUOTE_ALL,
            quotechar='"',
            escapechar='\\'
        )
        chunk_files.append(str(chunk_file))

        try:
            test_df = pd.read_csv(chunk_file, sep=",", quotechar='"', escapechar='\\', engine="python")
            if "smiles" not in test_df.columns:
                logger.error(f"{chunk_file} malformed: columns={test_df.columns.tolist()}")
        except Exception as e:
            logger.error(f"Failed to read {chunk_file}: {e}")

    logger.info(f"Generated {len(chunk_files)} chunk files")

    chunk_wrapper = ChunkWrapper(task_list, config)

    runner = ParallelWorkflowRunner(
        workflow_func=chunk_wrapper,
        config={"smiles_chunks": chunk_files},
        input_key="smiles_chunks",
        output_key="chunk_file",
        filename_pattern="{chunk_file}_processed.csv",
        combined_filename=config.get("output", {}).get("filename", "final_flexible_parallel_output.csv"),
        output_dir=str(chunk_dir),
        use_multiprocessing=True,
        reserve_cpus=4
    )

    chunk_results = runner.run()

    all_dfs = []
    all_mols = []

    if isinstance(chunk_results, pd.DataFrame):
        combined_df = chunk_results
    elif isinstance(chunk_results, (list, tuple)):
        for idx, result in enumerate(chunk_results):
            if not isinstance(result, tuple) or len(result) != 3:
                logger.warning(f"Invalid result format at index {idx}: {result}")
                continue
            chunk_id, df, mols = result
            if df is None or not isinstance(df, pd.DataFrame) or df.empty:
                continue
            all_dfs.append(df)
            if mols:
                all_mols.extend(mols)
        combined_df = pd.concat(all_dfs, ignore_index=True) if all_dfs else pd.DataFrame()
    else:
        combined_df = pd.DataFrame()

    logger.info(f"Combined dataframe shape: {combined_df.shape}")
    logger.debug(f"Combined dataframe columns: {combined_df.columns.tolist()}")

    final_csv = output_dir / config.get("output", {}).get("filename", "final_output.csv")

    try:
        combined_df.to_csv(
            final_csv,
            index=False,
            quoting=csv.QUOTE_ALL,
            quotechar='"',
            doublequote=True,
            escapechar=None,
            lineterminator='\n'
        )
        logger.info(f"Final CSV written: {final_csv}")
    except Exception as e:
        logger.error(f"Failed to write final CSV: {e}")

    if "smiles" in combined_df.columns:
        smi_file = final_csv.with_suffix(".smi")
        try:
            with open(smi_file, "w") as f:
                for smi in combined_df["smiles"].dropna():
                    f.write(smi + "\n")
            logger.info(f"SMILES file written: {smi_file}")
        except Exception as e:
            logger.error(f"Failed to write SMILES file: {e}")

    if all_mols:
        sdf_file = final_csv.with_suffix(".sdf")
        try:
            writer = Chem.SDWriter(str(sdf_file))
            for mol in all_mols:
                writer.write(mol)
            writer.close()
            logger.info(f"[dynamic_parallel_task_runner] SDF file written: {sdf_file}")
        except Exception as e:
            logger.error(f"Failed to write SDF file: {e}")
    else:
        logger.debug("No RDKit mols to write SDF.")

    if config.get("cleanup", True):
        try:
            shutil.rmtree(chunk_dir, ignore_errors=True)
        except Exception as e:
            logger.warning(f"Failed to cleanup chunk directory {chunk_dir}: {e}")

    return {"df": combined_df}

@register_workflow(
    "streamed_feature_runner",
    description="Stream large datasets (e.g. parquet), generate features per chunk, combine results, and run postprocessing analyses."
)
def streamed_feature_runner(config):
    input_file = config.get("input_file")
    if not input_file:
        raise ValueError("Missing 'input_file' in config")

    # === Load preprocess task list ===
    preprocess_tasks = config.get("workflow_preprocess", [])
    if not preprocess_tasks:
        raise ValueError("No 'workflow_preprocess' task list specified in config")

    # === Prepare output directories ===
    # Default chunk directory if nothing else is specified
    default_output_dir = Path("outputs/streamed_features")
    output_dir = default_output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    feature_files = []
    total_rows = 0
    max_row_groups = config.get("row_groups", None)

    # === Stream parquet row groups ===
    for i, (rg_index, df_chunk) in enumerate(
        tqdm(stream_parquet_rowgroups(input_file), desc="Row groups")
    ):
        if max_row_groups is not None and i >= max_row_groups:
            logger.info(f"[streamed_feature_runner] Stopping after {max_row_groups} row group(s) (prototype mode).")
            break

        logger.info(f"[streamed_feature_runner] Processing row group {rg_index}")
        if df_chunk is None or df_chunk.empty:
            logger.warning(f"Row group {rg_index} is empty, skipping.")
            continue

        current_data = {"df": df_chunk}

        # === Run all preprocessing tasks sequentially ===
        for task_name in preprocess_tasks:
            task_func = get_task(task_name)
            if task_func is None:
                registered = list_tasks()
                msg = (
                    f"[streamed_feature_runner] Task '{task_name}' not found.\n"
                    f"  → Did you mean one of these?\n"
                    f"    {registered}\n"
                    f"  Check your YAML 'workflow_preprocess:' section or @register_task names."
                )
                logger.error(msg)
                raise ValueError(msg)

            logger.info(f"[streamed_feature_runner] Running preprocess task '{task_name}'...")
            result = task_func(config, data=current_data)
            current_data = result if isinstance(result, dict) else {"df": result}

        df_features = current_data.get("df")
        if df_features is None or df_features.empty:
            logger.warning(f"Row group {rg_index} produced no features, skipping.")
            continue

        total_rows += len(df_features)

        # === Save per-chunk feature file ===
        # Store in generic streamed_features dir (fast local cache)
        chunk_dir = output_dir / "chunks"
        chunk_dir.mkdir(parents=True, exist_ok=True)

        feature_file = chunk_dir / f"features_chunk_{rg_index}.feather"
        df_features.to_feather(feature_file)
        feature_files.append(feature_file)

        del df_chunk, current_data, df_features
        gc.collect()

    logger.info(f"[streamed_feature_runner] Processed {total_rows:,} rows total across {len(feature_files)} chunks.")

    # === Combine chunked feather files into a single Parquet ===
    if feature_files:
        logger.info(f"[streamed_feature_runner] Combining {len(feature_files)} feature chunks into one Parquet file...")

        combined_df = pd.concat(
            [pd.read_feather(f) for f in feature_files],
            ignore_index=True
        )

        # Pull output config from the *last preprocessing task* (e.g. incremental_pca)
        last_task_name = preprocess_tasks[-1]
        last_task_cfg = config.get(last_task_name, {})
        output_cfg = last_task_cfg.get("output", {})

        combined_output_dir = Path(output_cfg.get("directory", "outputs/streamed_combined"))
        combined_output_dir.mkdir(parents=True, exist_ok=True)
        combined_output_filename = output_cfg.get("filename", "combined_output.parquet")
        combined_output_path = combined_output_dir / combined_output_filename

        combined_df.to_parquet(combined_output_path, index=False)
        logger.info(f"[streamed_feature_runner] Combined output written to: {combined_output_path}")
    else:
        logger.warning("[streamed_feature_runner] No feature chunks were generated; skipping combination step.")
        combined_output_path = None

    # === Run postprocessing tasks if defined ===
    post_tasks = config.get("workflow_postprocess", [])
    if post_tasks:
        logger.info(f"[streamed_feature_runner] Running postprocessing tasks: {post_tasks}")

        # Update config so downstream tasks see correct input file
        config = update_config_input_file(config, combined_output_path)

        current_data = {
            "input_file": str(combined_output_path),
            "feature_files": feature_files,
        }

        for task_name in post_tasks:
            task_func = get_task(task_name)
            if task_func is None:
                registered = list_tasks()
                msg = (
                    f"[streamed_feature_runner] Postprocessing task '{task_name}' not found.\n"
                    f"  → Did you mean one of these?\n"
                    f"    {registered}"
                )
                logger.error(msg)
                raise ValueError(msg)

            logger.info(f"[streamed_feature_runner] Running postprocess task '{task_name}'...")
            result = task_func(config, data=current_data)
            current_data = result if isinstance(result, dict) else {"df": result}

        logger.info("[streamed_feature_runner] All postprocessing tasks completed successfully.")
        return current_data

    # === Return summary if no postprocessing ===
    return {
        "feature_files": feature_files,
        "combined_output": str(combined_output_path) if combined_output_path else None
    }





