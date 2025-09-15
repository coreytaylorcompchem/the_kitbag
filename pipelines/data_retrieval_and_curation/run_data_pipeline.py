import argparse
import yaml
import os
from workflows import get_workflow
from pipeline.import_utilities import import_modules_from_dir

from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

def load_yaml(path):
    with open(path, 'r') as f:
        return yaml.safe_load(f)

def main():
    this_dir = os.path.dirname(os.path.abspath(__file__))
    import_modules_from_dir(os.path.join(this_dir, "workflows"), "workflows")
    import_modules_from_dir(os.path.join(this_dir, "modules"), "modules")

    # Load config from file
    parser = argparse.ArgumentParser()
    parser.add_argument('--params', required=True, help='Path to YAML config file')
    args = parser.parse_args()
    config = load_yaml(args.params)

    # Actually run the workflow
    workflow_name = config.get('workflow_name')
    if not workflow_name:
        raise ValueError("Missing 'workflow_name' in config")

    logger.info(f"Workflow name in config: '{workflow_name}'")

    workflow_func = get_workflow(workflow_name)
    if not workflow_func:
        raise ValueError(f"Workflow '{workflow_name}' not found")

    result_df = workflow_func(config)

    workflows_that_handle_output = {"chembl_multi_target"}

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
    # else:
    #     logger.info(f"Skipping write: Workflow '{workflow_name}' handles output internally.")

if __name__ == "__main__":
    main()
