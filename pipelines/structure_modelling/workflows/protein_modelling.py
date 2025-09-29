from workflows import register_workflow
from pipeline.logger import setup_logger
from pathlib import Path
import yaml
from backends import get_backend # later
from pipeline.task_registry import get_task

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_workflow("protein_modelling", description="Fix, model, and refine protein structure.")
def run(config_path: dict):
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    output_dir = Path(config.get("output_dir", "output"))
    output_dir.mkdir(parents=True, exist_ok=True)
    config["output_dir"] = output_dir

    backend_name = config["backend"]["name"]
    backend_kwargs = {k: v for k, v in config["backend"].items() if k != "name"}
    backend = get_backend(backend_name, **backend_kwargs)

    # Run tasks sequentially
    for step in config.get("workflow", []):
        logger.info(f"Running task: {step}")
        task_func = get_task(step)
        if not task_func:
            raise ValueError(f"Unknown task: {step}")
        task_func(backend, None, config)

    logger.info("Protein modeling pipeline complete.")
