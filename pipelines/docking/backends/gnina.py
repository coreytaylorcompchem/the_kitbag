import subprocess
from pathlib import Path
from typing import Tuple, Optional

from backends.base import BaseBackend
from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

class GninaBackend(BaseBackend):
    # manual list of supported tasks. Annoying but necessary for automatic discovery code.
    supported_tasks = [ 
        "adme_prediction",
        "standardise_ligand",
        "generate_conformers",
        "cluster_conformers",
        "save_final_conformers",
        "convert_to_pdbqt", 
        "prepare_receptor_pdbqt",
        "dock", 
        "active_learn_docking",
        "induced_fit_docking"]
    def __init__(self, executable_path: str = "gnina", use_gpu: bool = True):
        """
        Initialize Gnina backend.

        Parameters:
        -----------
        executable_path : str
            Path to gnina binary.
        use_gpu : bool
            If False, adds --cpu flag to force CPU mode.
        """
        super().__init__(executable_path=executable_path, use_gpu=use_gpu)

    def dock(self, ligand: dict, config: dict, **kwargs):
        """
        Run docking using GNINA for a single ligand.

        Expects:
        - receptor_pdbqt: in backend.cache["receptor_pdbqt"]
        - ligand_pdbqt: in ligand["pdbqt_path"]
        - docking center and size in config['docking']['center'] and ['size']
        - output_dir from config['output_dir']
        - Optional 'output_path' can be passed via kwargs
        """

        # Required paths
        receptor_path = self.cache.get("receptor_pdbqt")
        if receptor_path is None:
            raise ValueError("Receptor PDBQT not set in backend.cache['receptor_pdbqt'].")

        ligand_path = ligand.get("pdbqt_path")
        if ligand_path is None:
            raise ValueError(f"Missing 'pdbqt_path' for ligand '{ligand['name']}'.")

        # Use provided output_path or default
        output_path = kwargs.get("output_path")
        if output_path is None:
            output_path = Path(config["output_dir"]) / f"{ligand['name']}_docked.sdf"
        else:
            output_path = Path(output_path)

        # Docking box
        docking_cfg = config["docking"]
        center: Tuple[float, float, float] = docking_cfg["center"]
        size: Tuple[float, float, float] = docking_cfg["size"]
        num_modes = docking_cfg["n_output_binding_modes"]
        exhaustiveness = docking_cfg["exhaustiveness"]

        # Build GNINA command
        cmd = [
            str(self.executable_path),
            "-r", str(receptor_path),
            "-l", str(ligand_path),
            "--center_x", str(center[0]),
            "--center_y", str(center[1]),
            "--center_z", str(center[2]),
            "--size_x", str(size[0]),
            "--size_y", str(size[1]),
            "--size_z", str(size[2]),
            "-o", str(output_path),
            "--num_modes", str(num_modes),
            "--exhaustiveness", str(exhaustiveness)
        ]

        if not self.use_gpu:
            cmd.append("--cpu")

        extra_args: Optional[list] = docking_cfg.get("extra_args")
        if extra_args:
            cmd.extend(extra_args)

        # Run docking
        logger.debug(f"Running GNINA for ligand: {ligand['name']}")
        logger.debug(f"Command: {' '.join(cmd)}")

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            logger.error(f"GNINA failed:\n{result.stderr}")
            raise RuntimeError(f"GNINA docking failed for {ligand['name']}")

        logger.info(f"Docking completed for {ligand['name']}. Output: {output_path}")
        return [output_path]

