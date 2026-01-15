import sys
from pathlib import Path
import traceback
import json

import pandas as pd
from rdkit import Chem
import prolif as plf
from prolif import sdf_supplier

from workflows.utils.prolif import (
    select_best_poses,
    build_prolif_dataframe,
    plot_prolif_barcode,
)

from pipeline.logger import setup_logger

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


logger = setup_logger(
    "run_prolif",
    simple_format=True,
)

def main():
    """
    Standalone ProLIF runner.
    MUST be executed in a fresh Python process.
    """
    if len(sys.argv) not in [5, 6]:
        logger.error(
            "Usage: python -m pipeline.prolif_runner "
            "<docking_results.csv> <protein.pdb> <output_dir> [min_frequency]"
        )
        sys.exit(1)

    csv_path = Path(sys.argv[1])
    protein_pdb = Path(sys.argv[2])
    output_dir = Path(sys.argv[3])
    color_scheme = sys.argv[4] 
    interaction_colors = separated_interaction_colors if color_scheme == "separated" else grouped_interaction_colors
    min_freq = float(sys.argv[5]) if len(sys.argv) == 6 else 0.1

    output_dir.mkdir(exist_ok=True)

    try:
        best_df = select_best_poses(csv_path)

        fp_df = build_prolif_dataframe(
            best_df,
            protein_pdb=protein_pdb,
        )

        fp_df.to_csv(output_dir / "prolif_fingerprint.csv")

        plot_prolif_barcode(
            fp_df,
            output_dir,
            interaction_colors=interaction_colors, 
            min_frequency=min_freq,
        )

        logger.debug("ProLIF runner completed successfully.")

    except Exception as e:
        logger.error("ProLIF failed.")
        logger.error(str(e))
        traceback.print_exc()
        sys.exit(1) 


if __name__ == "__main__":
    main()
