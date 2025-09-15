from workflows import register_workflow
from pipeline.parallel_runner import ParallelWorkflowRunner
from pipeline.task_registry import get_task
from pipeline.logger import setup_logger
import pandas as pd
from pathlib import Path

logger = setup_logger(__name__)

def physchem_workflow_single(chunk_file: str):
    task = get_task("physchem_filtering")
    if not task:
        raise ValueError("Task 'physchem_filtering' not found")

    config = {"input_file": chunk_file}
    result = task(config)

    df = result.get("df") if result else None
    if df is None:
        df = pd.DataFrame()

    # Instead of returning dict, return tuple (filename, df)
    return (Path(chunk_file).stem, df)


@register_workflow(
    "physchem_filtering",
    description="Physchem filtering using calculated descriptors"
)
def physchem_filtering_parallel_workflow(config):
    input_file = config.get("input_file")
    if input_file is None:
        raise ValueError("config must supply 'input_file'")
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    df = pd.read_csv(input_file)
    if "smiles" not in df.columns:
        raise ValueError("CSV must contain a 'smiles' column")

    smiles_list = df["smiles"].dropna().unique().tolist()
    if not smiles_list:
        raise ValueError("No SMILES found in input file")

    chunk_size = config.get("chunk_size", 1000)
    output_base = Path(config.get("output", {}).get("directory", "outputs/physchem_parallel"))
    tmp_dir = output_base / "smiles_chunks"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    chunk_files = []
    for i in range(0, len(smiles_list), chunk_size):
        chunk = smiles_list[i:i+chunk_size]
        chunk_file = tmp_dir / f"smiles_chunk_{i // chunk_size}.csv"
        pd.DataFrame({"smiles": chunk}).to_csv(chunk_file, index=False)
        chunk_files.append(str(chunk_file))

    config_for_runner = dict(config)
    # Instead of "smiles_chunks" key, align with input_key param name for runner
    config_for_runner["smiles_chunks"] = chunk_files

    runner = ParallelWorkflowRunner(
        workflow_func=physchem_workflow_single,
        config={"smiles_chunks": chunk_files},
        input_key="smiles_chunks",
        output_key="chunk_file",  # this must match the format pattern below
        filename_pattern="{chunk_file}_filtered_df.csv",  # must match the tuple's first return value
        combined_filename=config.get("output", {}).get("filename", "filtered_physchem_combined.csv"),
        output_dir=str(output_base),
        use_multiprocessing=True,
        reserve_cpus=4
    )

    combined_result = runner.run()

    # Return combined dataframe as dict for consistency
    return {"df": combined_result}
