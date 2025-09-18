import re
import pandas as pd
import csv
import shutil
from pathlib import Path
from rdkit import Chem
from pipeline.task_registry import get_task
from pipeline.parallel_runner import ParallelWorkflowRunner
from pipeline.logger import setup_logger
from workflows import register_workflow

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

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

        logger.debug(f"[parallel_chunk_runner] Running {task_name} on chunk {chunk_file}")
        result = task(task_config, data=current_data)

        # Normalize result to a dict
        if isinstance(result, tuple):
            logger.debug(f"[parallel_chunk_runner] Task '{task_name}' returned tuple of length {len(result)}")

            if len(result) == 3:
                _, df_result, mols_result = result
                result = {"df": df_result, "mols": mols_result}
            elif len(result) == 2:
                _, df_result = result
                result = {"df": df_result, "mols": []}
            else:
                logger.debug(f"[parallel_chunk_runner] Task '{task_name}' returned unexpected tuple length {len(result)}")
                return (Path(chunk_file).stem, pd.DataFrame(), [])

        elif isinstance(result, dict):
            if "df" not in result or result["df"] is None or result["df"].empty:
                logger.debug(f"[parallel_chunk_runner] Task '{task_name}' returned empty or invalid dataframe")
                return (Path(chunk_file).stem, pd.DataFrame(), [])
        else:
            logger.error(f"[parallel_chunk_runner] Task '{task_name}' returned unsupported result type: {type(result)}")
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

    return (Path(chunk_file).stem, df, mols_out)

@register_workflow("dynamic_task_runner", description="Run a task list in parallel on input chunks and output CSV/SDF/SMI")
def dynamic_task_runner(config):
    input_file = config.get("input_file")
    if not input_file:
        raise ValueError("Missing 'input_file' in config")
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    task_list = config.get("workflow", [])
    if not task_list:
        raise ValueError("No 'workflow' list of tasks specified in config")

    df = pd.read_csv(input_file, sep=",", quotechar='"', escapechar='\\', engine="python")
    if "smiles" not in df.columns:
        raise ValueError("Input file must contain 'smiles' column")

    df["smiles"] = clean_smiles(df["smiles"])
    df = df[df["smiles"].str.len() > 0]

    smiles_list = df["smiles"].dropna().unique().tolist()
    chunk_size = config.get("chunk_size", 1000)

    output_dir = Path(config.get("output", {}).get("directory", "outputs/flexible_parallel"))
    output_dir.mkdir(parents=True, exist_ok=True)

    chunk_dir = output_dir / "chunks"
    if not chunk_dir.exists():
        chunk_dir.mkdir(parents=True, exist_ok=True)

    chunk_files = []
    for i in range(0, len(smiles_list), chunk_size):
        chunk = smiles_list[i: i + chunk_size]
        chunk_file = chunk_dir / f"chunk_{i//chunk_size}.csv"
        chunk_df = pd.DataFrame({"smiles": clean_smiles(pd.Series(chunk))})
        chunk_df.to_csv(chunk_file, index=False, quoting=csv.QUOTE_ALL, quotechar='"', escapechar='\\')
        chunk_files.append(str(chunk_file))

        # Optional check
        try:
            test_df = pd.read_csv(chunk_file, sep=",", quotechar='"', escapechar='\\', engine="python")
            if "smiles" not in test_df.columns or test_df.shape[1] != 1:
                logger.error(f"[Chunk Write Check] {chunk_file} malformed: columns={test_df.columns.tolist()}")
        except Exception as e:
            logger.error(f"[Chunk Write Check] Failed to read {chunk_file}: {e}")

    logger.info(f"[dynamic_parallel_task_runner] Generated {len(chunk_files)} chunk files")

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

    logger.info(f"[dynamic_parallel_task_runner] Combined dataframe shape: {combined_df.shape}")
    logger.info(f"[dynamic_parallel_task_runner] Combined dataframe columns: {combined_df.columns.tolist()}")

    final_csv = output_dir / config.get("output", {}).get("filename", "final_output.csv")

    try:
        combined_df.to_csv(
            final_csv,
            index=False,
            quoting=csv.QUOTE_ALL,
            quotechar='"',
            doublequote=True,  # ensure quotes inside fields are doubled
            escapechar=None,
            lineterminator='\n'
        )
        logger.info(f"[dynamic_parallel_task_runner] Final CSV written: {final_csv}")
    except Exception as e:
        logger.error(f"Failed to write final CSV: {e}")

    if "smiles" in combined_df.columns:
        smi_file = final_csv.with_suffix(".smi")
        try:
            with open(smi_file, "w") as f:
                for smi in combined_df["smiles"].dropna():
                    f.write(smi + "\n")
            logger.info(f"[dynamic_parallel_task_runner] SMILES file written: {smi_file}")
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
        logger.warning("[dynamic_parallel_task_runner] No RDKit mols to write SDF.")

    if config.get("cleanup", True):
        try:
            shutil.rmtree(chunk_dir, ignore_errors=True)
        except Exception as e:
            logger.warning(f"Failed to cleanup chunk directory {chunk_dir}: {e}")

    return {"df": combined_df}
