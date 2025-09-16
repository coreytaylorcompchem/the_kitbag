import argparse
import yaml
import os
import sys

from pipeline.import_utilities import import_modules_recursively
from workflows import get_workflow, load_all_workflows, list_workflows
from pipeline.logger import setup_logger

from pipeline.dependency_checker import check_dependencies, fail_if_missing

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def load_yaml(path):
    with open(path, 'r') as f:
        return yaml.safe_load(f)

def main():
    project_root = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, project_root)

    # Define what dependencies your pipeline expects (based on config)
    required_deps = {
        "openmm": "pip install openmm",
        "rdkit": "conda install -c conda-forge rdkit",
        "openff.toolkit": "conda install -c conda-forge openff-toolkit",
        "pdbfixer": "pip install pdbfixer",
        "tqdm": "pip install tqdm",
    }

    # Check and exit early if anything's missing
    missing = check_dependencies(required_deps)
    fail_if_missing(missing)

    # Import all modules and workflows
    import_modules_recursively(os.path.join(project_root, "modules"), "modules")
    import_modules_recursively(os.path.join(project_root, "workflows"), "workflows")
    load_all_workflows()

    logger.info(f"Registered workflows: {list_workflows()}")

    parser = argparse.ArgumentParser()
    parser.add_argument('--params', required=True, help='Path to YAML config file')
    args = parser.parse_args()

    config = load_yaml(args.params)

    platform = config["simulation"]["platform"]

    logger.info(f"Running on platform: {platform}")

    workflow_name = config.get("workflow_name")
    if not workflow_name:
        raise ValueError("Missing 'workflow_name' in config")

    logger.info(f"Running workflow: '{workflow_name}'")

    workflow_func = get_workflow(workflow_name)
    if not workflow_func:
        raise ValueError(f"Workflow '{workflow_name}' not found")

    result = workflow_func(config)

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
            logger.info(f"File exists and overwrite is false: {filepath}")
    else:
        logger.info(f"Workflow '{workflow_name}' handled its own output.")

if __name__ == "__main__":
    main()
