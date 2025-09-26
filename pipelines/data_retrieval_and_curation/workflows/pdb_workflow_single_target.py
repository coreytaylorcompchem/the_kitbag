from workflows import register_workflow
from pipeline.task_registry import get_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_workflow("pdb_single_target", description="Process PDBs for a single UniProt target.")
def run_pdb_single_target_workflow(config):
    steps = config.get("workflow", [])
    data = None

    for step in steps:
        task_func = get_task(step)
        if not task_func:
            raise ValueError(f"Task '{step}' not found in task registry.")
        logger.info(f">>>>>>>>>> Running task: {step}")
        data = task_func(config, data)

    logger.info("PDB workflow for single target complete.")
    return data
