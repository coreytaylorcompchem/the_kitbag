import pandas as pd
from pathlib import Path

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

from modules.utils.fpocket import run_fpocket_on_gpcr_structures

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_task("run_fpocket", category="Pocket detection", description="Run fpocket on protein structures to detect binding pockets.")
def run_fpocket(config: dict, data: dict = None) -> dict:
    """
    The main task to run fpocket on protein structures.
    """
    params = config.get("run_fpocket", {})

    gpcr_pdb_dir = Path(params.get("structure_directory"))
    output_dir = Path(params.get("output_directory"))
    max_pockets_per_structure = params.get("max_pockets_per_structure", 100)
    n_jobs = params.get("n_jobs", 1)
    
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Running fpocket on GPCR structures in {gpcr_pdb_dir}...")
    
    pocket_files = run_fpocket_on_gpcr_structures(
        gpcr_pdb_dir,
        output_dir,
        max_pockets_per_structure,
        n_jobs=n_jobs
    )
    
    if not pocket_files:
        logger.error("No pocket files were generated.")
        return {}

    # Combine all pocket CSVs into one
    all_pockets_df = pd.concat([pd.read_csv(f) for f in pocket_files], ignore_index=True)
    
    # Save combined pocket CSV
    combined_pocket_file = output_dir / "pockets.csv"
    all_pockets_df.to_csv(combined_pocket_file, index=False)
    logger.info(f"Combined pocket data saved to {combined_pocket_file}")

    return {"pockets_file": str(combined_pocket_file)}
