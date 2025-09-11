import os
import pandas as pd
from copy import deepcopy
from multiprocessing import Pool
from pipeline.task_registry import get_task
from pipeline.workflow_registry import register_workflow

def run_pipeline_for_target(args):
    """Run the configured workflow for a single UniProt ID. Must be top-level for multiprocessing."""
    uniprot_id, config = args
    output_dir = config.get("output", {}).get("directory", "outputs/parallel")
    os.makedirs(output_dir, exist_ok=True)

    local_config = deepcopy(config)
    local_config["uniprot_id"] = uniprot_id

    try:
        data = None
        for step in config["workflow"]:
            task_func = get_task(step)
            if not task_func:
                raise ValueError(f"Task '{step}' not found in registry.")
            data = task_func(local_config, data)

        if isinstance(data, dict) and "df" in data and isinstance(data["df"], pd.DataFrame):
            df = data["df"]
            if not df.empty:
                filename = f"{uniprot_id}_bioactivity.csv"
                path = os.path.join(output_dir, filename)
                df.to_csv(path, index=False)
                print(f"[{uniprot_id}] ✅ Saved: {path}")
                return {"df": df}
            else:
                print(f"[{uniprot_id}] ⚠️ No data after cleaning.")
        else:
            print(f"[{uniprot_id}] ⚠️ Invalid or empty data returned.")
    except Exception as e:
        print(f"❌ [{uniprot_id}] Error: {e}")

    # Always return a dict, even on failure
    return {"df": pd.DataFrame()}


def extract_valid_dfs(results):
    """Helper to extract non-empty DataFrames from list of task results."""
    return [
        r["df"] for r in results
        if isinstance(r, dict) and "df" in r and isinstance(r["df"], pd.DataFrame) and not r["df"].empty
    ]


@register_workflow(
    "chembl_multi_target",
    description="Retrieve and clean ChEMBL bioactivities for multiple targets."
)
def run_chembl_multi_target_parallel_workflow(config):
    uniprot_ids = config.get("uniprot_ids", [])
    if not uniprot_ids:
        raise ValueError("No UniProt IDs provided in config['uniprot_ids'].")

    output_config = config.get("output", {})
    output_dir = output_config.get("directory", "outputs/parallel")
    os.makedirs(output_dir, exist_ok=True)

    # Use N - 4 CPUs, minimum 1
    available_cpus = max(1, os.cpu_count() - 4)
    num_workers = min(available_cpus, len(uniprot_ids))
    print(f"🔧 Using {num_workers} workers to process {len(uniprot_ids)} targets")

    # Pair each UniProt ID with the full config
    args_list = [(uniprot_id, config) for uniprot_id in uniprot_ids]

    with Pool(processes=num_workers) as pool:
        results = pool.map(run_pipeline_for_target, args_list)

    valid_dfs = extract_valid_dfs(results)

    if not valid_dfs:
        print("❌ No valid dataframes to concatenate.")
        return pd.DataFrame()

    all_results = pd.concat(valid_dfs, ignore_index=True)

    combined_filename = output_config.get("filename", "combined_bioactivity.csv")
    combined_file = os.path.join(output_dir, combined_filename)
    all_results.to_csv(combined_file, index=False)

    print(f"\n✅ Combined data saved to {combined_file}")

    return all_results