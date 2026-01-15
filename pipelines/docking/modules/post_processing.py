from pathlib import Path
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

from workflows.utils.docking import (
    plot_docking_scatter_with_errorbars,
)

from workflows.utils.adme import (
    plot_adme_descriptor_grids,
)

from workflows.utils.prolif import (
    select_best_poses,
    build_prolif_dataframe,
    plot_prolif_barcode,
)

logger = setup_logger(__name__, simple_format=True)

separated_interaction_colors = { 
    "Hydrophobic": "#59e382", 
    "VdWContact": "#dfab43", 
    "HBAcceptor": "#59bee3", 
    "HBDonor": "#239fcd", 
    "XBAcceptor": "#ff9f02", 
    "XBDonor": "#ce8000", 
    "Cationic": "#e35959", 
    "Anionic": "#5979e3", 
    "CationPi": "#e359d8", 
    "PiCation": "#ea85e2", 
    "PiStacking": "#b559e3", 
    "EdgeToFace": "#c885ea", 
    "FaceToFace": "#a22ddc", 
    "MetalAcceptor": "#7da982", 
    "MetalDonor": "#609267", 
    "WaterBridge": "#323aa8", 
    }

grouped_interaction_colors = { 
    "Hydrophobic": "#59e382", 
    "VdWContact": "#dfab43", 
    "HBAcceptor": "#59bee3", 
    "HBDonor": "#239fcd", 
    "XBAcceptor": "#ff9f02", 
    "XBDonor": "#ce8000", 
    "Cationic": "#e35959", 
    "Anionic": "#5979e3", 
    "CationPi": "#e359d8", 
    "PiCation": "#ea85e2", 
    "PiStacking": "#b559e3", 
    "EdgeToFace": "#c885ea", 
    "FaceToFace": "#a22ddc", 
    "MetalAcceptor": "#7da982", 
    "MetalDonor": "#609267", 
    "WaterBridge": "#323aa8", 
    }

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


    # ProLIF step

    prolif_cfg = steps.get("prolif", {})
    if prolif_cfg.get("enabled", False):
        logger.info("Detecting P-L interactions and generating barcode plots.")

        min_freq = prolif_cfg.get("min_frequency", 0.1)
        color_scheme = prolif_cfg.get("color_scheme", "separated")

        interaction_colors = (
            separated_interaction_colors
            if color_scheme == "separated"
            else grouped_interaction_colors
        )

        best_df = select_best_poses(csv_path)

        try:
            fp_df = build_prolif_dataframe(
                best_df,
                protein_pdb=protein_pdb,
            )
        except Exception as e:
            logger.error(
                f"[post_processing] ProLIF failed unexpectedly: {e}"
            )
            fp_df = None

        if fp_df is None or fp_df.empty:
            logger.warning(
                "[post_processing] ProLIF skipped (no valid data)."
            )
        else:
            if prolif_cfg.get("output", {}).get("csv", False):
                fp_df.to_csv(analysis_dir / "prolif_fingerprint.csv")

            if prolif_cfg.get("output", {}).get("barcode", True):
                plot_prolif_barcode(
                    fp_df,
                    analysis_dir,
                    interaction_colors=interaction_colors,
                    min_frequency=min_freq,
                )

    # adme and docking score plots

    score_cfg = steps.get("docking_score_plots", {})
    if score_cfg.get("enabled", False):
        logger.info("Running docking score plots.")

        plot_docking_scatter_with_errorbars(csv_path, analysis_dir)

    adme_cfg = steps.get("adme_plots", {})
    if adme_cfg.get("enabled", False):
        logger.info("Running adme plots.")

        plot_adme_descriptor_grids(csv_path, analysis_dir)

    logger.info("Post_processing Completed.")

