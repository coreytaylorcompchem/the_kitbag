import subprocess
from pathlib import Path
from typing import Tuple
from backends.base import BaseBackend
from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

class UniDockBackend(BaseBackend):
    supported_tasks = [
        "standardise_ligand",
        "generate_conformers",
        "cluster_conformers",
        "save_final_conformers",
        "convert_to_pdbqt", 
        "prepare_receptor_pdbqt",
        "dock", ]
    def __init__(self, executable_path: str = "unidock", use_gpu: bool = True):
        super().__init__(executable_path=executable_path, use_gpu=use_gpu)

    def dock(self, ligand: dict, config: dict, **kwargs):
        """
        Dock ligand(s) using Uni-Dock.
        Supports both 'ensemble' and 'lowest_energy' modes.
        """

        receptor_path = self.cache.get("receptor_pdbqt")
        if receptor_path is None:
            raise ValueError("[Uni-Dock] receptor_pdbqt not found in backend.cache")

        docking_cfg = config["docking"]
        center: Tuple[float, float, float] = docking_cfg["center"]
        size: Tuple[float, float, float] = docking_cfg["size"]
        mode: str = docking_cfg.get("docking_mode", "ensemble")

        conformer_paths = ligand.get("pdbqt_paths", [])
        if not conformer_paths:
            raise ValueError(f"[Uni-Dock] No PDBQT paths found for ligand {ligand['name']}")

        # Handle docking mode
        if mode == "lowest_energy":
            conformer_paths = [conformer_paths[0]]

        successful_outputs = []

        for i, conf_path in enumerate(conformer_paths):
            conf_path = Path(conf_path)
            output_path = Path(config["output_dir"]) / f"{ligand['name']}_conf{i}_docked.pdbqt"

            cmd = [
                str(self.executable_path),
                "--receptor", str(receptor_path),
                "--ligand", str(conf_path),
                "--out", str(output_path),
                "--center_x", str(center[0]),
                "--center_y", str(center[1]),
                "--center_z", str(center[2]),
                "--size_x", str(size[0]),
                "--size_y", str(size[1]),
                "--size_z", str(size[2]),
                "--exhaustiveness", "8",
                "--num_modes", "20"
            ]

            if not self.use_gpu:
                logger.warning("[Uni-Dock] 'use_gpu' is set to False, but Uni-Dock does not support CPU-only mode directly. Ignoring.")

            logger.info("[Uni-Dock] Docking %s conf %d:", ligand['name'], i)
            logger.info("Running command: %s", " ".join(cmd))

            result = subprocess.run(cmd, capture_output=True, text=True)

            if result.returncode != 0:
                logger.error("[Uni-Dock] Failed docking %s conf %d:\n%s", ligand['name'], i, result.stderr)
                continue  # Optionally skip instead of raise

            successful_outputs.append(str(output_path))

        if not successful_outputs:
            raise RuntimeError(f"[Uni-Dock] Docking failed for all conformers of ligand {ligand['name']}")

        ligand["docked_outputs"] = successful_outputs
        logger.info("[Uni-Dock] Completed docking for %s (%d successful conformers)", ligand['name'], len(successful_outputs))
        return successful_outputs