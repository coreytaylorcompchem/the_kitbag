import os
import pandas as pd
from copy import deepcopy
from pipeline.task_registry import get_task
from pipeline.workflow_registry import register_workflow
from pipeline.parallel_runner import ParallelWorkflowRunner
import traceback
from tqdm import tqdm

def run_pipeline_for_target(local_config):
    uniprot_id = local_config.get("uniprot_id")
    output_dir = local_config.get("output", {}).get("directory", "outputs/tox_targets")
    os.makedirs(output_dir, exist_ok=True)

    try:
        data = None
        workflow_steps = local_config.get("workflow", [])
        for step in workflow_steps:
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

        # compile final df
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
        print(traceback.format_exc())
        return pd.DataFrame()


@register_workflow("chembl_tox_targets", description="Retrieve bioactivity data for tox-relevant targets")
def run_chembl_tox_targets_parallel_workflow(config):
    runner = ParallelWorkflowRunner(
        workflow_func=run_pipeline_for_target,
        config=config,
        input_key="uniprot_ids",       # list of uniprot_ids in config
        output_key="uniprot_id",       # each run gets a single uniprot_id in config['uniprot_id']
        output_dir=config.get("output", {}).get("directory", "outputs/tox_targets"),
        filename_pattern="{uniprot_id}_bioactivity.csv",
        combined_filename=config.get("output", {}).get("filename", "combined_tox_bioactivity.csv"),
        use_multiprocessing=True,
        reserve_cpus=4
    )
    
    return runner.run()
