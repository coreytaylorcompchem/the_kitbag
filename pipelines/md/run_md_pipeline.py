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

    required_deps = {
        "openmm": "pip install openmm",
        "rdkit": "conda install -c conda-forge rdkit",
        "openff.toolkit": "conda install -c conda-forge openff-toolkit",
        "pdbfixer": "conda install -c conda-forge pdbfixer",
        "tqdm": "pip install tqdm",
    }

    missing = check_dependencies(required_deps)
    fail_if_missing(missing)

    import_modules_recursively(os.path.join(project_root, "modules"), "modules")
    import_modules_recursively(os.path.join(project_root, "workflows"), "workflows")
    load_all_workflows()

    logger.info(f"Registered workflows: {list_workflows()}")

    parser = argparse.ArgumentParser()
    parser.add_argument('--params', required=True, help='Path to YAML config file')
    args = parser.parse_args()

    config = load_yaml(args.params)

    platform = config.get("setup_simulation", {}).get("platform", "CPU")
    logger.info(f"Running on platform: {platform}")

    workflow_name = config.get("workflow_name")
    if not workflow_name:
        raise ValueError("Missing 'workflow_name' in config")

    logger.info(f"Running workflow: '{workflow_name}'")
    workflow_func = get_workflow(workflow_name)
    if not workflow_func:
        raise ValueError(f"Workflow '{workflow_name}' not found")

    # Hack: inject flags to control workflow steps based on YAML list
    steps = config.get("workflow", [])
    config["run_prepare_system"] = "prepare_system" in steps
    config["run_setup_simulation"] = "setup_simulation" in steps
    config["run_minimize"] = "minimize" in steps
    config["run_heat_and_equilibrate"] = "heat_and_equilibrate" in steps
    config["run_production"] = "production" in steps

    workflow_result = workflow_func(config)

    if isinstance(workflow_result, dict) and "df" in workflow_result:
        result_df = workflow_result["df"]
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
