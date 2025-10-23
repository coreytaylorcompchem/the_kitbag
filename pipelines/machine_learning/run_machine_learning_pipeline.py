import argparse
import os
import sys
import yaml
from pathlib import Path

from workflows import get_workflow, load_all_workflows, list_workflows
from pipeline.import_utilities import import_modules_recursively
from pipeline.task_registry import list_tasks, get_task

from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

def parse_args():
    parser = argparse.ArgumentParser(description="Run docking pipeline via dynamically selected workflow.")
    parser.add_argument(
        "-p", "--params", required=True,
        help="Path to the YAML configuration file."
    )
    parser.add_argument(
        "--debug", action="store_true",
        help="Enable debug output (optional)."
    )
    return parser.parse_args()


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

    args = parse_args()
    config_path = Path(args.params)

    if not config_path.exists():
        logger.error(f"Config file not found: {config_path}")
        sys.exit(1)

    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    workflow_name = config.get("workflow_name")
    if not workflow_name:
        logger.error("Config must include a 'workflow_name' key.")
        sys.exit(1)

    workflow_func = get_workflow(workflow_name)
    if workflow_func is None:
        logger.error(f"Unknown workflow: '{workflow_name}'")
        logger.info("Available workflows:")
        for name in list_workflows():
            logger.info(f"  - {name}")
        sys.exit(1)

    if args.debug:
        logger.debug(f"Selected workflow: {workflow_name}")
        logger.debug(f"Config path: {config_path.resolve()}")

    try:
        workflow_func(str(config_path))
    except Exception as e:
        logger.error(f"[FATAL] Workflow '{workflow_name}' failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
