from workflows import register_workflow
from modules.md import MDWorkflow
from modules.md_post_processing import MDPostProcessingWorkflow
from backends.openmm_backend import OpenMMBackend
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_workflow(
    "molecular_dynamics",
    description="Perform molecular dynamics simulation using specified backend."
)
def run_basic_md_workflow(config: dict):
    logger.info(f"Running workflow: {config.get('workflow_name', 'unknown')}")
    logger.info(f"Workflow steps: {config.get('workflow', [])}")

    backend = OpenMMBackend(config)
    workflow = MDWorkflow(backend, config)
    post_workflow = MDPostProcessingWorkflow(config)

    # Step mapping
    step_mapping = {
        "prepare_system": workflow.prepare_system,
        "setup_simulation": workflow.setup_simulation,
        "minimize": workflow.minimize,
        "heat_and_equilibrate": workflow.heat_and_equilibrate,
        "production": workflow.production,
        "post_processing": post_workflow.post_processing,
    }

    context = {}

    for task_name in config.get("workflow", []):
        task_func = step_mapping.get(task_name)
        if not task_func:
            logger.warning(f"Task '{task_name}' not found in MDWorkflow or MDPostProcessingWorkflow")
            continue

        logger.info(f"Running task '{task_name}'...")
        result = task_func()
        context[task_name] = result
        logger.info(f"Task '{task_name}' completed.")

    return context
