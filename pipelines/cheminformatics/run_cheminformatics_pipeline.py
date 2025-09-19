import argparse
import yaml
import os
import sys
import pandas as pd

from pipeline.import_utilities import import_modules_recursively
from pipeline.task_registry import list_tasks
from pipeline.logger import setup_logger
from workflows import get_workflow, load_all_workflows, list_workflows

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

def load_yaml(path):
    with open(path, 'r') as f:
        return yaml.safe_load(f)

def main():
    project_root = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, project_root)

    # # Import all workflows recursively
    # import_modules_recursively(
    #     base_dir=os.path.join(project_root, "workflows"),
    #     base_package="workflows"
    # )
    load_all_workflows()  # if this does extra workflow discovery
    logger.info(f"Discovered workflows: {list_workflows()}")

    # Import all tasks recursively
    import_modules_recursively(
        base_dir=os.path.join(project_root, "modules"),
        base_package="modules"
    )

    logger.info(f"Registered tasks: {list_tasks()}")

    # Parse CLI args
    parser = argparse.ArgumentParser()
    parser.add_argument('--params', required=True, help='Path to YAML config file')
    args = parser.parse_args()

    # Load YAML config
    config = load_yaml(args.params)

    # Retrieve workflow name from config
    workflow_name = config.get("workflow_name")
    if not workflow_name:
        raise ValueError("Missing 'workflow_name' in config")

    logger.info(f"Running workflow: '{workflow_name}'")

    # Get the registered workflow function
    workflow_func = get_workflow(workflow_name)
    if not workflow_func:
        raise ValueError(f"Workflow '{workflow_name}' not found")

    # Run the workflow
    result = workflow_func(config)

    # Check for valid result structure
    if isinstance(result, dict):
        if "df" in result:
            result_df = result["df"]
        else:
            logger.error(f"Workflow result does not contain 'df' key: {result}")
            result_df = pd.DataFrame()  # fallback to empty DataFrame
    else:
        logger.error(f"Unexpected result type: {type(result)}")
        result_df = pd.DataFrame()  # fallback to empty DataFrame

    # Output handling for workflows that do NOT write their own output
    workflows_that_handle_output = {
        "physchem_filtering"  # these manage their own file writing
    }

    if workflow_name not in workflows_that_handle_output:
        if result_df is not None and hasattr(result_df, 'to_csv'):
            output_cfg = config.get("output", {})
            filename = output_cfg.get("filename", "output.csv")
            directory = output_cfg.get("directory", "outputs")
            os.makedirs(directory, exist_ok=True)
            filepath = os.path.join(directory, filename)

            if not os.path.exists(filepath) or output_cfg.get("overwrite", False):
                result_df.to_csv(filepath, index=False)
                logger.info(f"\nOutput saved to {filepath}")
            else:
                logger.info(f"\nOutput file exists and overwrite is false: {filepath}")
    else:
        logger.info(f"Workflow '{workflow_name}' handles output internally.")

if __name__ == "__main__":
    main()
