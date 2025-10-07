import argparse
import yaml
import os
import sys
import pandas as pd

from pipeline.import_utilities import import_modules_recursively
from pipeline.task_registry import list_tasks, get_task
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

    load_all_workflows()
    logger.info(f"Discovered workflows: {list_workflows()}")

    import_modules_recursively(
        base_dir=os.path.join(project_root, "modules"),
        base_package="modules"
    )
    logger.info(f"Registered tasks: {list_tasks()}")

    parser = argparse.ArgumentParser()
    parser.add_argument('--params', required=True, help='Path to YAML config file')
    args = parser.parse_args()

    config = load_yaml(args.params)

    workflow_name = config.get("workflow_name")
    if not workflow_name:
        raise ValueError("Missing 'workflow_name' in config")

    logger.info(f"Running workflow: '{workflow_name}'")

    workflow_func = get_workflow(workflow_name)
    if not workflow_func:
        raise ValueError(f"Workflow '{workflow_name}' not found")

    # Handle dynamic_task_runner with parallel and serial workflows
    if workflow_name == "dynamic_task_runner":
        # Run parallel workflow tasks if specified
        if "workflow_parallel" in config and config["workflow_parallel"]:
            config_parallel = dict(config)
            config_parallel["workflow"] = config_parallel["workflow_parallel"]
            logger.info(f"Running parallel workflow tasks: {config_parallel['workflow']}")
            result_parallel = workflow_func(config_parallel)

            # Prepare config for serial workflow tasks
            parallel_output_dir = config_parallel.get("output", {}).get("directory", "outputs/flexible_parallel")
            parallel_output_file = config_parallel.get("output", {}).get("filename", "final_output.csv")
            combined_output_path = os.path.join(parallel_output_dir, parallel_output_file)

            config_serial = dict(config)
            config_serial["workflow"] = config_serial.get("workflow_serial", [])
            config_serial["chunk_size"] = 0  # disable chunking for serial workflow
            config_serial["input_file"] = combined_output_path

            if config_serial["workflow"]:
                logger.info(f"Running serial workflow tasks: {config_serial['workflow']}")
                result_serial = workflow_func(config_serial)
                result = result_serial
            else:
                result = result_parallel
        else:
            # No parallel workflow: run serial workflow only
            config_single = dict(config)
            config_single["workflow"] = config_single.get("workflow_serial", [])
            if not config_single["workflow"]:
                raise ValueError("No 'workflow_parallel' or 'workflow_serial' specified in config")
            logger.info(f"Running workflow tasks: {config_single['workflow']}")
            result = workflow_func(config_single)
    else:
        # Other workflows: run as usual
        result = workflow_func(config)

    # Handle output for workflows that don't manage their own output
    workflows_that_handle_output = {
        "physchem_filtering"  # example, adjust as needed
    }

    if workflow_name not in workflows_that_handle_output:
        if isinstance(result, dict) and "df" in result:
            result_df = result["df"]
            output_cfg = config.get("output", {})
            filename = output_cfg.get("filename", "output.csv")
            directory = output_cfg.get("directory", "outputs")
            os.makedirs(directory, exist_ok=True)
            filepath = os.path.join(directory, filename)

            if not os.path.exists(filepath) or output_cfg.get("overwrite", False):
                result_df.to_csv(filepath, index=False)
                logger.info(f"Output saved to {filepath}")
            else:
                logger.info(f"Output file exists and overwrite is False: {filepath}")
        else:
            logger.warning("Workflow result does not contain a valid 'df' key")
    else:
        logger.info(f"Workflow '{workflow_name}' handles output internally.")

if __name__ == "__main__":
    main()
