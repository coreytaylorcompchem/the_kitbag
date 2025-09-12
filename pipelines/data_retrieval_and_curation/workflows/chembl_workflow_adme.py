import os
import pandas as pd
from copy import deepcopy
from pipeline.task_registry import get_task
from pipeline.workflow_registry import register_workflow
from pipeline.parallel_runner import ParallelWorkflowRunner

from workflows.utils import process_readout_dataframe


def run_adme_pipeline_for_readout(config: dict) -> dict:
    """
    Run pipeline for a single ADME readout specified in config['readout'].

    Args:
        config (dict): Config dict containing at least the 'readout' key.

    Returns:
        dict: {'readout': <readout_name>, 'df': <resulting_dataframe>}
    """
    readout = str(config.get("readout"))
    if not readout:
        raise ValueError("No readout specified in config['readout'].")

    local_config = deepcopy(config)

    readout_filters = local_config.get("filters", {}).get(readout)
    if readout_filters:
        local_config.update(readout_filters)

    data = None
    for step in local_config.get("workflow", []):
        task_func = get_task(step)
        if not task_func:
            raise ValueError(f"Task '{step}' not found in registry.")
        data = task_func(local_config, data)

    if isinstance(data, dict) and "df" in data:
        df = data["df"]
    else:
        df = data

    # Apply check for empty dataframes from querying here

    if isinstance(data, dict) and "df" in data:
        df = data["df"]
    else:
        df = data

    df = process_readout_dataframe(readout, df)

    return {"readout": readout, "df": df}

# Checks the YAML file for mistakes and errors out if there's missing readout data retrieve (e.g. no assay_type)

def validate_filters_section(config):
    filters = config.get("filters", {})
    readouts = config.get("readout", [])
    if isinstance(readouts, str):
        readouts = [readouts]

    for readout, filter_set in filters.items():
        if readout not in readouts:
            print(f"[Info] Skipping validation for non-readout filter key: '{readout}'")
            continue

        if not isinstance(filter_set, dict):
            raise ValueError(f"[Config Error] Filter for '{readout}' must be a dictionary, got {type(filter_set)}.")

        if "assay_type" not in filter_set:
            raise ValueError(f"[Config Error] Filter for readout '{readout}' is missing required 'assay_type'.")


@register_workflow("chembl_adme_data", description="Retrieve ADME data from ChEMBL (parallelized by readout)")
def run_chembl_adme_workflow(config):
    readouts = config.get("readout", [])
    if isinstance(readouts, str):
        readouts = [readouts]
    if not readouts:
        raise ValueError("No readouts specified in config['readout'].")
    
    validate_filters_section(config)  # Here's where the YAML check is run

    local_config = deepcopy(config)
    local_config["readout"] = readouts

    runner = ParallelWorkflowRunner(
        workflow_func=run_adme_pipeline_for_readout,
        config=local_config,
        input_key="readout",
        output_key="readout",
        output_dir=local_config.get("output", {}).get("directory", "outputs/adme"),
        filename_pattern="{readout}_adme.csv",
        combined_filename=local_config.get("output", {}).get("filename", "combined_adme_data.csv"),
    )

    return runner.run()
