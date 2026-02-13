import os

from pathlib import Path
import subprocess
from backends.base import BaseBackend
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

class GlideBackend(BaseBackend):
    supported_tasks = [
        "standardise_ligand",
        "protonate_ligand",
        "prepare_ligand_glide",
        "prepare_receptor_grid",
        "dock",
        "compute_rdkit_descriptors",
        "adme_prediction",
    ]

    def __init__(self, executable_path: str, **kwargs):
        if executable_path is None:
            raise ValueError("Glide backend requires 'executable_path' in backend config.")
        
        super().__init__(executable_path=executable_path, **kwargs)
        self.schrodinger_home = Path(executable_path).resolve()
        if not self.schrodinger_home.exists():
            raise EnvironmentError(f"Schrödinger path does not exist: {self.schrodinger_home}")

        self.schrodinger_home = Path(executable_path).resolve()
        if not self.schrodinger_home.exists():
            raise EnvironmentError(f"Schrödinger path does not exist: {self.schrodinger_home}")

    def _exe(self, name):
        # default: top-level executables
        exe = os.path.join(self.executable_path, name)
        
        if name.lower() == "prepwizard":
            exe = os.path.join(self.executable_path, "utilities", "prepwizard")
        
        if not os.path.exists(exe):
            raise FileNotFoundError(f"{name} not found at {exe}")
        return exe

    def dock(self, ligand: dict, config: dict, output_path=None, **kwargs):
        """Dock a single ligand (may have multiple conformers). Always returns list of Path/str."""
        grid_file = self.cache.get("receptor_grid")
        ligand_input = ligand.get("glide_input")
        if not grid_file or not ligand_input:
            raise ValueError("Receptor grid or ligand input missing for Glide docking.")

        output_dir = Path(config["output_dir"])
        output_dir.mkdir(parents=True, exist_ok=True)

        n_poses = config.get("docking", {}).get("n_output_binding_modes", 10)

        # Ensure ligand_input is list
        ligand_inputs = ligand_input if isinstance(ligand_input, (list, tuple)) else [ligand_input]

        docked_outputs = []

        for idx, ligand_file in enumerate(ligand_inputs):
            # Generate unique output path per conformer
            if output_path and len(ligand_inputs) == 1:
                output_file = Path(output_path)
            elif output_path:
                output_file = Path(output_path).with_name(
                    f"{Path(output_path).stem}_conf{idx}{Path(output_path).suffix}"
                )
            else:
                output_file = output_dir / f"{ligand['name']}_conf{idx}_docked.maegz"

            cmd = [
                self._exe("glide"),
                "-HOST", "localhost",
                "-IN", str(ligand_file),
                "-GRID", str(grid_file),
                "-OUT", str(output_file),
                "-POSES", str(n_poses),
            ]

            logger.info(f"[GlideBackend] Docking {ligand['name']} (conf {idx})...")
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                logger.error(f"Glide docking failed:\n{result.stderr}")
                continue

            # Convert MAEGZ → SDF
            output_sdf = output_file.with_suffix(".sdf")
            convert_cmd = [self._exe("structconvert"), str(output_file), str(output_sdf)]
            convert_result = subprocess.run(convert_cmd, capture_output=True, text=True)
            if convert_result.returncode != 0:
                logger.error(f"Structconvert failed:\n{convert_result.stderr}")
                continue

            docked_outputs.append(output_sdf)

        if not docked_outputs:
            raise RuntimeError(f"Glide docking failed for {ligand['name']}")

        # Store string paths for downstream tasks
        ligand["docked_outputs"] = [str(p) for p in docked_outputs]
        return {"docked_output": [str(p) for p in docked_outputs]}
