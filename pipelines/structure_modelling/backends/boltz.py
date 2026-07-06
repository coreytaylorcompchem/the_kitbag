import os
import yaml
import subprocess
import logging
import json
from pathlib import Path
from backends.base import BaseStructureTool

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def validate_msa_policy(
    sequences: dict,
    msas: dict,
    use_msa_server: bool,
    tool_name: str,
):
    missing = [
        chain_id
        for chain_id in sequences
        if not msas or chain_id not in msas
    ]

    if missing and not use_msa_server:
        raise ValueError(
            f"{tool_name}: use_msa_server=false, but no local MSA was provided "
            f"for chains: {', '.join(missing)}"
        )

def generate_boltz_yaml(
    sequences: dict,
    ligands: list,
    templates: list,
    msas: dict,
    yaml_path: Path,
):

    data = {
        "sequences": []
    }

    # templates structures, if present.
    
    if templates:
        data["templates"] = []

        for t in templates:
            # detect extension
            ext = Path(t).suffix.lower()

            if ext == ".cif":
                entry = {"cif": t}
            elif ext == ".pdb":
                entry = {"pdb": t}
            else:
                raise ValueError(f"Unsupported template format: {t}")

            data["templates"].append(entry)
               
        logger.info(
            f"Templates: {len(templates) if templates else 0}"
        )

    # proteins 
    for chain_id, seq in sequences.items():
        protein_entry = {
            "id": chain_id,
            "sequence": seq,
        }

        if msas and chain_id in msas:
            protein_entry["msa"] = str(msas[chain_id])

        data["sequences"].append({
            "protein": protein_entry
        })

    # ligands 
    for i, smiles in enumerate(ligands):
        lig_id = f"L{i+1}"

        data["sequences"].append({
            "ligand": {
                "id": lig_id,
                "smiles": smiles
            }
        })

    # strict validation
    if not data["sequences"]:
        raise ValueError("No valid sequences or ligands provided for Boltz input")

    yaml_str = yaml.dump(data, sort_keys=False)

    yaml_path.write_text(yaml_str)

def log_parameters(params: dict, title="Parameters"):
    if not params:
        logger.info(f"[{title}] No parameters")
        return

    key_width = max(len(str(k)) for k in params.keys()) + 2

    logger.info(f"[{title}]")
    logger.info("─" * (key_width + 20))

    for k, v in params.items():
        logger.info(f"{k:<{key_width}} {v}")

    logger.info("─" * (key_width + 20))

class BoltzBackend(BaseStructureTool):
    name = "boltz"
    
    def __init__(self, config=None):
            super().__init__()
            self.config = config or {}
    def run(
        self,
        run_id: int,
        device: int,
        output_dir: Path,
        sequences: dict,
        ligands: list = None,
        templates: list = None,
        msas: dict = None
    ):
        
        use_msa_server = bool(inf_cfg.get("use_msa_server", False))

        validate_msa_policy(
            sequences=sequences,
            msas=msas or {},
            use_msa_server=use_msa_server,
            tool_name=self.name,
        )

        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = str(device)

        run_dir = output_dir / f"run_{run_id}"
        run_dir.mkdir(parents=True, exist_ok=True)

        yaml_file = run_dir / "input.yaml"
     
        generate_boltz_yaml(
            sequences,
            ligands or [],
            templates or [],
            msas or {},
            yaml_file
        )

        out_dir = run_dir / "results"
        
        cmd = [
            "boltz",
            "predict",
            str(yaml_file),
            "--out_dir", str(out_dir),
        ]
        
        # config extraction
        sp_cfg = self.config.get("structure_prediction", {})
        inf_cfg = sp_cfg.get("inference", {})

        # accelerator 
        accelerator = inf_cfg.get("accelerator", "gpu")

        # Add yaml parameters to Boltz command.
        # accelerator
        cmd += ["--accelerator", accelerator]

        # optional flags
        if inf_cfg.get("use_msa_server", True):
            cmd.append("--use_msa_server")

        if inf_cfg.get("no_kernels", True):
            cmd.append("--no_kernels")            
        if inf_cfg.get("use_potentials", False):
            cmd.append("--use_potentials")

        # numeric params
        if "diffusion_samples" in inf_cfg:
            cmd += ["--diffusion_samples", str(inf_cfg["diffusion_samples"])]

        if "sampling_steps" in inf_cfg:
            cmd += ["--sampling_steps", str(inf_cfg["sampling_steps"])]

        if "step_scale" in inf_cfg:
            cmd += ["--step_scale", str(inf_cfg["step_scale"])]

        if "recycling_steps" in inf_cfg:
            cmd += ["--recycling_steps", str(inf_cfg["recycling_steps"])]

        if "sampling_steps_affinity" in inf_cfg:
            cmd += ["--sampling_steps_affinity", str(inf_cfg["sampling_steps_affinity"])]

        if "diffusion_samples_affinity" in inf_cfg:
            cmd += ["--diffusion_samples_affinity", str(inf_cfg["diffusion_samples_affinity"])]

        # clean run label
        run_label = f"Boltz run={run_id} device={device}"

        # runtime parameters
        runtime_params = {
            "accelerator": accelerator,
            "use_msa_server": inf_cfg.get("use_msa_server", True),
            "no_kernels": inf_cfg.get("no_kernels", True),
            "use_potentials": inf_cfg.get("use_potentials", True),
        }

        # sampling parameters
        sampling_params = {
            "diffusion_samples": inf_cfg.get("diffusion_samples"),
            "sampling_steps": inf_cfg.get("sampling_steps"),
            "step_scale": inf_cfg.get("step_scale"),
            "recycling_steps": inf_cfg.get("recycling_steps"),
        }

        # affinity parameters 
        affinity_params = {
            "diffusion_samples_affinity": inf_cfg.get("diffusion_samples_affinity"),
            "sampling_steps_affinity": inf_cfg.get("sampling_steps_affinity"),
        }

        # remove None values
        sampling_params = {k: v for k, v in sampling_params.items() if v is not None}
        affinity_params = {k: v for k, v in affinity_params.items() if v is not None}

        # log cleanly
        log_parameters(runtime_params, f"{run_label} | Runtime")
        log_parameters(sampling_params, f"{run_label} | Sampling")

        if affinity_params:
            log_parameters(affinity_params, f"{run_label} | Affinity")

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

        # find boltz_results directory
        boltz_dirs = list(out_dir.glob("boltz_results_*"))

        if not boltz_dirs:
            logger.warning(f"No boltz_results directory found in {out_dir}")
            return []

        samples = []

        # recursively find JSON files
        json_files = list(boltz_dirs[0].rglob("confidence*.json"))

        if not json_files:
            logger.warning(f"No confidence JSON files found under {boltz_dirs[0]}")
            return []

        for json_file in json_files:
            logger.debug(f"Found JSON: {json_file}")

            data = json.loads(json_file.read_text())

            # find corresponding structure file

            parent_dir = json_file.parent

            # confidence_input_model_3.json
            # -> input_model_3

            model_stem = json_file.stem.replace(
                "confidence_",
                ""
            )
    
            structure_file = None

            for ext in [".cif", ".pdb"]:

                candidate = parent_dir / f"{model_stem}{ext}"

                if candidate.exists():
                    structure_file = candidate
                    break

            if structure_file is None:
                logger.warning(
                    f"No structure file found for {json_file}"
                )
                continue

            samples.append({
                "structure": str(structure_file),
                "plddt": (
                    data.get("complex_plddt")
                    or data.get("mean_plddt")
                    or data.get("plddt")
                ),
                "iptm": data.get("iptm"),
                "ptm": data.get("ptm"),
                "iplddt": data.get("complex_iplddt"),
            })

        return samples

