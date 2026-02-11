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

    def dock(self, ligand: dict, config: dict, **kwargs):
        """
        Dock ligand(s) using Glide.
        Supports both single-ligand and conformer ensembles if ligand['glide_input'] is a list.
        Stores docked output path(s) in ligand['docked_outputs'].
        """
        grid_file = self.cache.get("receptor_grid")
        ligand_input = ligand.get("glide_input")
        if not grid_file or not ligand_input:
            raise ValueError("Receptor grid or ligand input missing for Glide docking.")

        output_dir = Path(config["output_dir"])
        output_dir.mkdir(parents=True, exist_ok=True)

        docking_mode = config.get("docking", {}).get("docking_mode", "ensemble").lower()
        n_poses = config.get("docking", {}).get("n_output_binding_modes", 10)

        # Allow glide_input to be either a single file or list (ensemble)
        if isinstance(ligand_input, (list, tuple)):
            ligand_inputs = ligand_input
        else:
            ligand_inputs = [ligand_input]

        docked_outputs = []

        for idx, ligand_file in enumerate(ligand_inputs):
            output_file = output_dir / f"{ligand['name']}_conf{idx}_docked.maegz"
            cmd = [
                "glide",
                "-HOST", "localhost",
                "-IN", str(ligand_file),
                "-GRID", str(grid_file),
                "-OUT", str(output_file),
                "-POSES", str(n_poses),
            ]
            logger.info(f"[GlideBackend] Docking {ligand['name']} (conf {idx}) using Glide...")
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                logger.error(f"Glide docking failed for {ligand['name']} (conf {idx}):\n{result.stderr}")
                continue
            docked_outputs.append(output_file)

        if not docked_outputs:
            raise RuntimeError(f"Glide docking failed for all conformers of {ligand['name']}")

        ligand["docked_outputs"] = [str(p) for p in docked_outputs]
        logger.info(f"[GlideBackend] Completed docking for {ligand['name']}, outputs: {docked_outputs}")
        
        return docked_outputs if len(docked_outputs) > 1 else docked_outputs[0]
