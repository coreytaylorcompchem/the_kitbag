import argparse
import yaml
import os
import sys

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

    # Import all workflows and tasks recursively
    load_all_workflows()
    logger.info(f"Discovered workflows: {list_workflows()}")

    import_modules_recursively(
        base_dir=os.path.join(project_root, "modules"),
        base_package="modules"
    )
    logger.info(f"Registered tasks: {list_tasks()}")

    # CLI argument parsing
    parser = argparse.ArgumentParser(description="Run a specified workflow with YAML config.")
    parser.add_argument('--params', required=True, help='Path to YAML config file')
    args = parser.parse_args()

    # Load config
    config = load_yaml(args.params)

    workflow_name = config.get("workflow_name")
    if not workflow_name:
        raise ValueError("Missing 'workflow_name' in config")

    logger.info(f"Running workflow: '{workflow_name}'")

    workflow_func = get_workflow(workflow_name)
    if not workflow_func:
        raise ValueError(f"Workflow '{workflow_name}' not found")

    # Run workflow with config dictionary
    result = workflow_func(config)

    # Basic handling for workflows that return dict results with 'df' keys
    if isinstance(result, dict) and "df" in result:
        df = result["df"]
        output_cfg = config.get("output", {})
        filename = output_cfg.get("filename", "output.csv")
        directory = output_cfg.get("directory", "outputs")
        os.makedirs(directory, exist_ok=True)
        filepath = os.path.join(directory, filename)

        if not os.path.exists(filepath) or output_cfg.get("overwrite", False):
            df.to_csv(filepath, index=False)
            logger.info(f"Output saved to {filepath}")
        else:
            logger.info(f"Output file exists and overwrite is false: {filepath}")
    else:
        logger.info("Workflow completed. No tabular output to save or output handled internally.")

if __name__ == "__main__":
    main()
