from workflows import register_workflow
from pipeline.parallel_runner import ParallelWorkflowRunner
from pipeline.task_registry import get_task
from pipeline.logger import setup_logger

import pandas as pd

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

def workflow_function(config):
    steps = config.get("workflow", [])
    data = None
    for step in steps:
        task_func = get_task(step)
        if not task_func:
            raise ValueError(f"Task '{step}' not found in registry.")
        data = task_func(config, data)
    return data

@register_workflow(
    "chembl_multi_target",
    description="Retrieve, standardise and collate bioactivities for multiple targets - CHEMBL."
)
def run_chembl_multi_target_parallel_workflow(config):
    steps = config.get("workflow", [])

    # Case 1: UniProt IDs already provided
    if config.get("uniprot_ids"):
        logger.info("✅ Using provided UniProt IDs from config.")

    # Case 2: Try to retrieve UniProt IDs using a task
    elif "retrieve_protein_class_target_list" in steps:
        logger.info("🔍 Attempting to retrieve UniProt IDs via 'retrieve_protein_class_target_list'...")

        task_func = get_task("retrieve_protein_class_target_list")
        if not task_func:
            raise ValueError("Task 'retrieve_protein_class_target_list' not found in registry.")

        task_func(config)  # updates config['uniprot_ids']

        if not config.get("uniprot_ids"):
            logger.error("🛑 No UniProt IDs found after running 'retrieve_protein_class_target_list'.")
            return pd.DataFrame()

        # Remove the task since it's already run
        steps = [step for step in steps if step != "retrieve_protein_class_target_list"]
        config["workflow"] = steps

    # Case 3: Neither provided nor retrievable
    else:
        raise ValueError("❌ No UniProt IDs found or retrievable. Provide 'uniprot_ids' or include 'retrieve_protein_class_target_list' in workflow.")

    # Proceed with parallel bioactivity fetching
    runner = ParallelWorkflowRunner(
        workflow_func=workflow_function,
        config=config,
        input_key="uniprot_ids",
        output_key="uniprot_id",
        filename_pattern="{uniprot_id}_bioactivity.csv",
        combined_filename=config.get("output", {}).get("filename", "combined_bioactivity.csv"),
        output_dir=config.get("output", {}).get("directory", "outputs/parallel")
    )

    return runner.run()
