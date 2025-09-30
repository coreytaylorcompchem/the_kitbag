from pipeline.task_registry import register_task
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

@register_task("build_peptide", category="Peptide modelling", description="Build peptides from sequence file")
def build_peptide(backend, config, **kwargs):
    logger.info("Building peptides...")
    results = backend.run()  # build + minimization done here
    return results

@register_task("minimise_peptide", category="Peptide modelling", description="Minimize peptide structures")
def minimise_peptide(backend, config, **kwargs):
    # Since minimization is done inside backend.run(), this can be a noop or just log
    logger.info("Peptides minimized (included in build_peptide step).")
    return {}
