import os
import sys
import argparse
import yaml
import importlib
import logging

from modules import load_all_tasks
from workflows import load_all_workflows
from pipeline.task_registry import get_task
from workflows import get_workflow
from pipeline.backend_registry import get_backend_class

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def main():
    parser = argparse.ArgumentParser(description="Quantum chemistry pipeline runner")
    parser.add_argument('--params', required=True, help="Path to YAML config file")
    args = parser.parse_args()

    # Load config YAML
    with open(args.params, 'r') as f:
        config = yaml.safe_load(f)

    # Load input geometry file
    input_file = config.get('input_file')
    if not input_file or not os.path.exists(input_file):
        raise FileNotFoundError("Input file not found or not specified in config (under 'input_file')")

    # Import all tasks and workflows
    load_all_tasks()
    load_all_workflows()

    task_type = config.get('task', '').lower()
    if not task_type:
        raise ValueError("Missing 'task' in config.")

    if task_type == "workflow":
        workflow_name = config.get('workflow_name')
        if not workflow_name:
            raise ValueError("Missing 'workflow_name' in config.")

        workflow_func = get_workflow(workflow_name)
        if not workflow_func:
            raise ValueError(f"No workflow named '{workflow_name}' found.")

        logger.info(f"Running workflow: {workflow_name}")
        workflow_func(input_file, config)  # pass xyz path and full config

    else:
        # Single-task execution (non-workflow)
        task_func = get_task(task_type)
        if not task_func:
            raise ValueError(f"Unknown task: {task_type}")
        
        task_config = config.get(task_type, config)  # fallback to global config
        backend_name = task_config.get('backend')
        if not backend_name:
            raise ValueError(f"Missing 'backend' key in config for task '{task_type}'")

        BackendClass = get_backend_class(backend_name)
        backend = BackendClass(config)

        logger.info(f"Running task '{task_type}' with backend '{backend_name}'")
        task_func(backend, input_file, task_config)


if __name__ == '__main__':
    main()
