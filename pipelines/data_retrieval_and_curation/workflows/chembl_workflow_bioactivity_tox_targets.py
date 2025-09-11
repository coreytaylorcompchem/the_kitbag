import os
import pandas as pd
from copy import deepcopy
from multiprocessing import Pool
from pipeline.task_registry import get_task
from pipeline.workflow_registry import register_workflow

import traceback
from tqdm import tqdm

def run_pipeline_for_target(args):
    uniprot_id, config = args
    local_config = deepcopy(config)
    local_config["uniprot_id"] = uniprot_id

    output_dir = config.get("output", {}).get("directory", "outputs/tox_targets")
    os.makedirs(output_dir, exist_ok=True)

    try:
        data = None
        # Optionally, progress over steps for a target
        workflow_steps = config.get("workflow", [])
        for i, step in enumerate(workflow_steps):
            task_func = get_task(step)
            if not task_func:
                raise ValueError(f"Task '{step}' not found in registry.")

            if step == "retrieve_chembl_bioactivities":
                result = task_func(local_config, data)
                data = result

            elif step == "clean_bioactivities":
                readouts = local_config.get("readout", ["IC50"])
                if isinstance(readouts, str):
                    readouts = [readouts]

                cleaned_frames = []
                for readout in tqdm(readouts, desc=f"[{uniprot_id}] Cleaning readouts", leave=False):
                    clean_config = deepcopy(local_config)
                    clean_config["readout"] = [readout]

                    input_df = data if isinstance(data, pd.DataFrame) else data.get("df")
                    cleaned_result = task_func(clean_config, input_df)
                    cleaned_df = cleaned_result.get("df") if isinstance(cleaned_result, dict) else cleaned_result

                    if cleaned_df is not None and not cleaned_df.empty:
                        cleaned_df = cleaned_df.copy()
                        cleaned_df["readout"] = readout.upper()
                        print(f"[{uniprot_id}] Cleaned data for '{readout}': shape = {cleaned_df.shape}")
                        cleaned_frames.append(cleaned_df)
                    else:
                        print(f"[{uniprot_id}] No data for readout '{readout}'")

                if cleaned_frames:
                    combined_cleaned = pd.concat(cleaned_frames, ignore_index=True)
                    print(f"[{uniprot_id}] Combined cleaned data shape: {combined_cleaned.shape}")
                    data = {"df": combined_cleaned, "readout": None}
                else:
                    print(f"[{uniprot_id}] No readout data cleaned for any requested readouts.")
                    data = {"df": pd.DataFrame(), "readout": None}

            elif step == "retrieve_compound_data":
                try:
                    result = task_func(local_config, data)
                    data = result
                except ValueError as e:
                    if "Bioactivity data is empty" in str(e) or "no readout selected" in str(e):
                        print(f"[{uniprot_id}] Skipping 'retrieve_compound_data': {e}")
                        data = {"df": pd.DataFrame(), "readout": None}
                    else:
                        raise

            else:
                result = task_func(local_config, data)
                data = result

        # After all steps, extract final DataFrame to save
        if isinstance(data, dict) and "df" in data:
            final_df = data["df"]
        elif isinstance(data, pd.DataFrame):
            final_df = data
        else:
            raise ValueError(f"Unexpected final output type: {type(data)}")

        if not final_df.empty:
            file_path = os.path.join(output_dir, f"{uniprot_id}_bioactivity.csv")
            final_df.to_csv(file_path, index=False)
            print(f"[{uniprot_id}] Saved: {file_path}")
            return final_df
        else:
            print(f"⚠️⚠️⚠️ [{uniprot_id}] There was no bioactivity data for {uniprot_id}")
            return pd.DataFrame()

    except Exception as e:
        print(f"❌ [{uniprot_id}] Error: {e}")
        print(traceback.format_exc())  # Print full traceback for easier debugging
        return pd.DataFrame()

@register_workflow("chembl_tox_targets", description="Retrieve bioactivity data for tox-relevant targets")
def run_chembl_tox_targets_parallel_workflow(config):
    uniprot_ids = config.get("uniprot_ids", [])
    output_cfg = config.get("output", {})
    output_dir = output_cfg.get("directory", "outputs/tox_targets")
    os.makedirs(output_dir, exist_ok=True)

    # 🔧 Reserve 4 CPUs
    available_cpus = max(1, os.cpu_count() - 4)
    num_workers = min(available_cpus, len(uniprot_ids))
    print(f"Using {num_workers} workers to process {len(uniprot_ids)} targets")

    args_list = [(uid, config) for uid in uniprot_ids]

    results = []
    with Pool(processes=num_workers) as pool:
        # Use imap_unordered for tqdm progress bar support
        for result in tqdm(pool.imap_unordered(run_pipeline_for_target, args_list), total=len(args_list), desc="Targets"):
            results.append(result)

    valid_results = [df for df in results if not df.empty]

    if not valid_results:
        print("❗ No valid bioactivity data was collected from any target.")
        return pd.DataFrame()

    all_data = pd.concat(valid_results, ignore_index=True)

    if not all_data.empty:
        combined_path = os.path.join(output_dir, output_cfg.get("filename", "combined_tox_bioactivity.csv"))
        all_data.to_csv(combined_path, index=False)
        print(f"\nCombined tox bioactivity data saved to {combined_path}")

    return all_data
