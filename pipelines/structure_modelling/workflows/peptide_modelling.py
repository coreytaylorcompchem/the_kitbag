from workflows import register_workflow
from pipeline.logger import setup_logger
from pathlib import Path
from backends.peptidebuilder import PeptideBuilderBackend
from pipeline.task_registry import get_task

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_workflow("peptide_modelling", 
                   description="Build, minimise peptides and calculate metrics.")
def run(config: dict):
    output_dir = Path(config.get("output_dir", "output"))
    output_dir.mkdir(parents=True, exist_ok=True)
    config["output_dir"] = output_dir

    backend_name = config["backend"]["name"]
    
    peptide_cfg = config.get("build_peptide", {})
    sequence_file = peptide_cfg.get("sequence_file")
    cap_termini = peptide_cfg.get("cap_termini", True)

    if backend_name == "peptidebuilder":
        backend = PeptideBuilderBackend(
            sequence_file=sequence_file,
            output_dir=output_dir,
            cap_termini=cap_termini,
        )
    else:
        raise ValueError(f"Unsupported backend: {backend_name}")

    context = {}

    for step in config.get("workflow", []):
        logger.info(f"Running task: {step}")
        task_func = get_task(step)
        if not task_func:
            raise ValueError(f"Unknown task: {step}")

        result = task_func(backend, config, **context)

        if step in {"generate_peptide_conformers", "optimise_peptide_conformers"}:
            context["energy_map"] = result
        elif step == "analyse_peptide_flexibility":
            context["flexibility"] = result

    logger.info("Peptide modelling workflow complete.")
