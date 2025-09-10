import os
import pandas as pd
from copy import deepcopy
from multiprocessing import Pool, cpu_count
from pipeline.task_registry import get_task
from pipeline.workflow_registry import register_workflow


def run_pipeline_for_target(args):
    """Must be at top level to be picklable by multiprocessing"""
    uniprot_id, config = args
    from_path = config.get("output", {}).get("directory", "outputs/parallel")
    os.makedirs(from_path, exist_ok=True)

    local_config = deepcopy(config)
    local_config["uniprot_id"] = uniprot_id

    try:
        data = None
        for step in config["workflow"]:
            task_func = get_task(step)
            if not task_func:
                raise ValueError(f"Task '{step}' not found in registry.")
            data = task_func(local_config, data)

        if data is not None and not data.empty:
            filename = f"{uniprot_id}_bioactivity.csv"
            path = os.path.join(from_path, filename)
            data.to_csv(path, index=False)
            print(f"[{uniprot_id}] Saved: {path}")
            return data
        else:
            print(f"[{uniprot_id}] No data for {uniprot_id}")
            return pd.DataFrame()

    except Exception as e:
        print(f"❌ [{uniprot_id}] Error: {e}")
        return pd.DataFrame()


@register_workflow(
    'chembl_multi_target',
    description="Retrieve and clean ChEMBL bioactivities for multiple targets."
)
def run_chembl_multi_target_parallel_workflow(config):
    uniprot_ids = config.get("uniprot_ids", [])
    output_config = config.get("output", {})
    output_dir = output_config.get("directory", "outputs/parallel")
    os.makedirs(output_dir, exist_ok=True)

    num_workers = min(cpu_count(), len(uniprot_ids))
    print(f"Using {num_workers} workers to process {len(uniprot_ids)} targets")

    # Pass both the target and the config to each worker
    args_list = [(uniprot_id, config) for uniprot_id in uniprot_ids]

    with Pool(processes=num_workers) as pool:
        results = pool.map(run_pipeline_for_target, args_list)

    all_results = pd.concat([df for df in results if not df.empty], ignore_index=True)

    combined_filename = output_config.get("filename", "combined_bioactivity.csv")
    combined_file = os.path.join(output_dir, combined_filename)

    all_results.to_csv(combined_file, index=False)
    
    print(f"\nCombined data saved to {combined_file}")

    return all_results