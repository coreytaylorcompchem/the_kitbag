import argparse
import yaml
import os
from pipeline.workflow_registry import get_workflow
from pipeline.import_utilities import import_modules_from_dir

def load_yaml(path):
    with open(path, 'r') as f:
        return yaml.safe_load(f)

def main():
    # Dynamically import all tasks and workflows before using registries
    this_dir = os.path.dirname(os.path.abspath(__file__))
    import_modules_from_dir(os.path.join(this_dir, "workflows"), "workflows")
    import_modules_from_dir(os.path.join(this_dir, "tasks"), "tasks")

    # Load config
    parser = argparse.ArgumentParser()
    parser.add_argument('--params', required=True, help='Path to YAML config file')
    args = parser.parse_args()
    config = load_yaml(args.params)

    # Get and run workflow
    workflow_name = config.get('workflow_name')
    if not workflow_name:
        raise ValueError("Missing 'workflow_name' in config")

    workflow_func = get_workflow(workflow_name)
    if not workflow_func:
        raise ValueError(f"Workflow '{workflow_name}' not found")

    result_df = workflow_func(config)

    # Handle output
    if result_df is not None and hasattr(result_df, 'to_csv'):
        output_cfg = config.get("output", {})
        filename = output_cfg.get("filename", "output.csv")
        directory = output_cfg.get("directory", "outputs")
        os.makedirs(directory, exist_ok=True)
        filepath = os.path.join(directory, filename)

        if not os.path.exists(filepath) or output_cfg.get("overwrite", False):
            result_df.to_csv(filepath, index=False)
            print(f"\n✅ Output saved to {filepath}")
        else:
            print(f"\n⚠️ Output file exists and overwrite is false: {filepath}")

if __name__ == "__main__":
    main()
