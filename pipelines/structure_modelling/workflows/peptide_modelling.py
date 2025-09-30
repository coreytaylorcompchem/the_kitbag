from workflows import register_workflow
from pipeline.logger import setup_logger
from pathlib import Path
from backends.peptidebuilder import PeptideBuilderBackend
from pipeline.task_registry import get_task

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_workflow("peptide_modelling", description="Build and minimise peptides.")
def run(config: dict):
    output_dir = Path(config.get("output_dir", "output"))
    output_dir.mkdir(parents=True, exist_ok=True)
    config["output_dir"] = output_dir

    backend_name = config["backend"]["name"]
    if backend_name != "build_peptide":
        raise ValueError(f"Unsupported backend: {backend_name}")

    minimization_cfg = config.get("minimization", {})
    peptide_cfg = config.get("peptide", {})
    sequence_file = peptide_cfg.get("sequence_file")

    backend = PeptideBuilderBackend(sequence_file=sequence_file, output_dir=output_dir, minimization_cfg=minimization_cfg)

    for step in config.get("workflow", []):
        logger.info(f"Running task: {step}")
        task_func = get_task(step)
        if not task_func:
            raise ValueError(f"Unknown task: {step}")
        task_func(backend, config)

    logger.info("Peptide modelling workflow complete.")
