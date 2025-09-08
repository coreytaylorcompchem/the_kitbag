import argparse
import sys
import yaml
from pathlib import Path

from workflows import get_workflow, list_workflows


def parse_args():
    parser = argparse.ArgumentParser(description="Run docking pipeline via dynamically selected workflow.")
    parser.add_argument(
        "-c", "--config", required=True,
        help="Path to the YAML configuration file."
    )
    parser.add_argument(
        "--debug", action="store_true",
        help="Enable debug output (optional)."
    )
    return parser.parse_args()


def main():
    args = parse_args()
    config_path = Path(args.config)

    if not config_path.exists():
        print(f"[ERROR] Config file not found: {config_path}")
        sys.exit(1)

    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    workflow_name = config.get("workflow")
    if not workflow_name:
        print("[ERROR] Config must include a 'workflow' key.")
        sys.exit(1)

    workflow_func = get_workflow(workflow_name)
    if workflow_func is None:
        print(f"[ERROR] Unknown workflow: '{workflow_name}'")
        print("Available workflows:")
        for name in list_workflows():
            print(f"  - {name}")
        sys.exit(1)

    if args.debug:
        print(f"[DEBUG] Selected workflow: {workflow_name}")
        print(f"[DEBUG] Config path: {config_path.resolve()}")

    try:
        workflow_func(str(config_path))
    except Exception as e:
        print(f"[FATAL] Workflow '{workflow_name}' failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
