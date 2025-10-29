import yaml
from pathlib import Path

from workflows import register_workflow
from modules.md import MDWorkflow
from backends.openmm_backend import OpenMMBackend
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_workflow(
    "molecular_dynamics",
    description="Perform molecular dynamics simulation using OpenMM."
)
def run_basic_md_workflow(config: dict):
    """
    Run MD workflow using tasks defined in MDWorkflow.
    Expects `config` to be a dict already loaded from YAML.
    """
    logger.info(f"Running workflow: {config.get('workflow_name', 'unknown')}")
    logger.info(f"Workflow steps: {config.get('workflow', [])}")

    # Initialize backend and workflow instance
    backend = OpenMMBackend(config)
    workflow = MDWorkflow(backend, config)

    # Map workflow step names to MDWorkflow instance methods
    step_mapping = {
        "prepare_system": workflow.prepare_system,
        "setup_simulation": workflow.setup_simulation,
        "minimize": workflow.minimize,
        "heat_and_equilibrate": workflow.heat_and_equilibrate,
        "production": workflow.production,
    }

    # Execute each task in the YAML-defined workflow
    for task_name in config.get("workflow", []):
        task_func = step_mapping.get(task_name)
        if not task_func:
            logger.warning(f"Task '{task_name}' not found in MDWorkflow")
            continue

        logger.info(f"Running task '{task_name}'...")
        task_func()
        logger.info(f"Task '{task_name}' completed.")

    return {}
