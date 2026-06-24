from workflows import register_workflow
from pipeline.task_registry import get_task
from pathlib import Path
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=True, simple_format=True)


@register_workflow(
    "structure_prediction",
    description="Parallel AI-based structure prediction"
)
def run(config, **kwargs):
    output_dir = Path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    class PipelineState:
        def __init__(self):
            self.cache = {}

    state = PipelineState()
    steps = config.get("workflow", [])

    if not steps:
        logger.error("No workflow steps defined in config. Please make sure to add tasks to the yaml.")
        return state

    for step in steps:
        logger.info(f"[Workflow] Running step: {step}")

        task = get_task(step)

        if task is None:
            raise ValueError(f"Task '{step}' not found in registry")

        task(state, config)

    logger.info("[Workflow] Completed all steps")

    return state