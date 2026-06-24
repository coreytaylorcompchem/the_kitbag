import os
import yaml
import subprocess
import logging
import json
from pathlib import Path
from backends.base import BaseStructureTool

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=True, simple_format=True)

def generate_boltz_yaml(sequences: dict, yaml_path: Path):
    """
    sequences = {"A": "SEQ1", "B": "SEQ2"}
    """

    data = {
        "version": 1,
        "sequences": [],
        "targets": [
            {
                "name": "complex",
                "sequences": list(sequences.keys())
            }
        ]
    }

    for chain_id, seq in sequences.items():
        data["sequences"].append({
            "protein": {
                "id": chain_id,
                "sequence": seq
            }
        })

    yaml_path.write_text(yaml.dump(data))

class BoltzBackend(BaseStructureTool):
    name = "boltz"

    def run(
        self,
        run_id: int,
        device: int,
        output_dir: Path,
        sequences: dict
    ):
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = str(device)

        run_dir = output_dir / f"run_{run_id}"
        run_dir.mkdir(parents=True, exist_ok=True)

        yaml_file = run_dir / "input.yaml"

        generate_boltz_yaml(sequences, yaml_file)

        out_dir = run_dir / "results"

        cmd = [
            "boltz",
            "predict",
            str(yaml_file),
            "--out_dir", str(out_dir),
            "--accelerator", "gpu",
            "--use_msa_server",
            "--diffusion_samples", "3",
            "--sampling_steps", "150",
            "--step_scale", "1.4",
            "--no_kernels",
        ]

        logger.info(f"[Boltz] run={run_id} gpu={device}")
        subprocess.run(cmd, check=True, env=env)

        results = self._parse_results(out_dir)
        
        if not results:
            logger.warning(f"[WARNING] No samples returned for run {run_id}")


        return {
            "run_id": run_id,
            "device": device,
            "results": results,
            "tool": self.name,
        }
    
    def _parse_results(self, out_dir: Path):

        # --- find boltz_results directory ---
        boltz_dirs = list(out_dir.glob("boltz_results_*"))

        if not boltz_dirs:
            logger.warning(f"No boltz_results directory found in {out_dir}")
            return []

        samples = []

        # --- recursively find JSON files ---
        json_files = list(boltz_dirs[0].rglob("confidence*.json"))

        if not json_files:
            logger.warning(f"No confidence JSON files found under {boltz_dirs[0]}")
            return []

        for json_file in json_files:
            logger.info(f"[DEBUG] Found JSON: {json_file}")

            data = json.loads(json_file.read_text())

            # --- find matching structure file in same dir ---
            parent_dir = json_file.parent

            structure_file = None
            for ext in ["*.cif", "*.pdb"]:
                files = list(parent_dir.glob(ext))
                if files:
                    structure_file = files[0]
                    break

            if structure_file is None:
                logger.warning(f"No structure file found next to {json_file}")
                continue

            samples.append({
                "structure": str(structure_file),
                "plddt": data.get("mean_plddt"),
                "ptm": data.get("ptm"),
                "iptm": data.get("iptm"),
            })

        return samples

