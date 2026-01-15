from pathlib import Path
import subprocess
import sys
import os
import json

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

from workflows.utils.docking import (
    plot_docking_scatter_with_errorbars,
)

from workflows.utils.adme import (
    plot_adme_descriptor_grids,
)

logger = setup_logger(__name__, simple_format=True)

@register_task(
    "post_processing",
    category="Post-processing",
    description="Run post-docking analysis steps (ProLIF, plots, etc.)",
)
def post_processing(backend, ligand, config, **kwargs):
    """
    Dataset-level task.
    Executed once after docking.
    """

    # -----------------------------------------
    # Skip per-ligand calls
    # -----------------------------------------
    if ligand is not None:
        return None

    cfg = config.get("post_processing", {})
    if not cfg.get("enabled", False):
        logger.info("[post_processing] Disabled.")
        return None

    output_dir = Path(config["output_dir"])
    analysis_dir = output_dir / "analysis"
    analysis_dir.mkdir(exist_ok=True)

    csv_path = output_dir / "docking_results.csv"
    protein_pdb = output_dir / "receptor_protonated.pdb"

    if not csv_path.exists():
        raise FileNotFoundError(csv_path)
    if not protein_pdb.exists():
        raise FileNotFoundError(protein_pdb)

    steps = cfg.get("steps", {})

    
    # ProLIF workaround (Prolif is not thread-safe otherwise)
    
    prolif_cfg = steps.get("prolif", {})
    min_freq = prolif_cfg.get("min_frequency", 0.1)
    color_scheme = prolif_cfg.get("color_scheme", "separated")

    if prolif_cfg.get("enabled", False):
        
        logger.info("Detecting P–L interactions.")

        cmd = [
            sys.executable,
            "-u",             
            "-m",
            "pipeline.prolif_runner",
            str(csv_path),
            str(protein_pdb),
            str(analysis_dir),
            color_scheme,
            str(min_freq),
        ]

        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = "1"
        env["MKL_NUM_THREADS"] = "1"

        project_root = Path(__file__).resolve().parents[1]
        env["PYTHONPATH"] = str(project_root)

        result = subprocess.run(
            cmd,
            # stdout=subprocess.PIPE,
            # stderr=subprocess.PIPE,
            # text=True,
            env=env,
        )

        if result.returncode != 0:
            logger.warning("ProLIF failed.")
            if result.stderr:
                logger.warning(result.stderr.strip())
        else:
            logger.info("ProLIF completed successfully.")


    # Docking score plots

    score_cfg = steps.get("docking_score_plots", {})
    if score_cfg.get("enabled", False):
        logger.info("Running docking score plots.")
        plot_docking_scatter_with_errorbars(csv_path, analysis_dir)


    # ADME plots

    adme_cfg = steps.get("adme_plots", {})
    if adme_cfg.get("enabled", False):
        logger.info("Running adme plots.")
        plot_adme_descriptor_grids(csv_path, analysis_dir)

    logger.info("Post_processing Completed.")
