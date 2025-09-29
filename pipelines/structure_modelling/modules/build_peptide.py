from pipeline.task_registry import register_task
from pathlib import Path
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_task("build_peptide_batch", category="Peptide modeling", description="Predict structures for multiple peptides using ColabFold")
def build_peptide_batch(backend, ligand, config, **kwargs):
    from pathlib import Path

    output_dir = Path(config["output_dir"])
    sequence_file = Path(config["peptide"]["sequence_file"])

    result = backend.build_peptides_from_file(sequence_file, output_dir)
    backend.cache["peptide_models"] = {name: str(pdb_path) for name, pdb_path in result}

    logger.info(f"Finished modeling {len(result)} peptides with ColabFold.")

