import os
import yaml
import subprocess
import json
import statistics

from pathlib import Path

from backends.base import BaseStructureTool
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)


def generate_intellifold_yaml(
    sequences: dict,
    ligands: list,
    templates: list,
    yaml_path: Path,
):
    """
    IntelliFold YAML schema is similar in spirit to Boltz:
      sequences:
        - protein:
            id: A
            sequence: ...
        - ligand:
            id: L1
            smiles: ...
    """

    data = {
        "sequences": []
    }

    # --- proteins ---
    for chain_id, seq in sequences.items():
        data["sequences"].append({
            "protein": {
                "id": chain_id,
                "sequence": seq,
            }
        })

    # --- ligands ---
    for i, smiles in enumerate(ligands or []):
        lig_id = f"L{i+1}"

        data["sequences"].append({
            "ligand": {
                "id": lig_id,
                "smiles": smiles,
            }
        })

    # --- templates ---
    # IntelliFold supports template usage via CLI --use_template.
    # For precomputed template paths, schema may differ slightly by version,
    # so keep this simple for now and validate with one generated YAML.
    if templates:
        data["templates"] = []

        for t in templates:
            ext = Path(t).suffix.lower()

            if ext == ".cif":
                data["templates"].append({"cif": str(t)})
            elif ext == ".pdb":
                data["templates"].append({"pdb": str(t)})
            else:
                raise ValueError(f"Unsupported template format: {t}")

    if not data["sequences"]:
        raise ValueError("No valid sequences or ligands provided for IntelliFold input")

    yaml_path.write_text(
        yaml.dump(data, sort_keys=False)
    )


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


def _mean_if_list(value):
    if isinstance(value, list) and value:
        return float(statistics.mean(value))
    return value


def _first_present(*values):
    for v in values:
        if v is not None:
            return v
    return None


class IntelliFoldBackend(BaseStructureTool):
    name = "intellifold"

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
    ):
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = str(device)

        run_dir = output_dir / f"intellifold_run_{run_id}"
        run_dir.mkdir(parents=True, exist_ok=True)

        yaml_file = run_dir / "input.yaml"

        out_dir = run_dir / "results"

        sp_cfg = self.config.get("structure_prediction", {})
        inf_cfg = sp_cfg.get("inference", {})

        use_template = bool(inf_cfg.get("use_template", False))

        available_templates = templates or []
        effective_templates = available_templates

        if available_templates and not use_template:
            logger.info(
                f"Templates were provided for run {run_id}, but use_template=false; "
                "IntelliFold template usage will be disabled."
            )
            effective_templates = []

        generate_intellifold_yaml(
            sequences=sequences,
            ligands=ligands or [],
            templates=effective_templates,
            yaml_path=yaml_file,
        )

        cmd = [
            "intellifold",
            "predict",
            str(yaml_file),
            "--out_dir",
            str(out_dir),
        ]

        # --- common IntelliFold options ---
        if "model" in inf_cfg:
            cmd += ["--model", str(inf_cfg["model"])]

        if "seed" in inf_cfg:
            cmd += ["--seed", str(inf_cfg["seed"])]

        # Accept both IntelliFold-native and generic pipeline names
        num_diffusion_samples = inf_cfg.get(
            "num_diffusion_samples",
            inf_cfg.get("diffusion_samples"),
        )

        recycling_iters = inf_cfg.get(
            "recycling_iters",
            inf_cfg.get("recycling_steps"),
        )

        sampling_steps = inf_cfg.get("sampling_steps")

        if num_diffusion_samples is not None:
            cmd += [
                "--num_diffusion_samples",
                str(num_diffusion_samples),
            ]

        if sampling_steps is not None:
            cmd += [
                "--sampling_steps",
                str(sampling_steps),
            ]

        if recycling_iters is not None:
            cmd += [
                "--recycling_iters",
                str(recycling_iters),
            ]

        if "precision" in inf_cfg:
            cmd += [
                "--precision",
                str(inf_cfg["precision"]),
            ]

        if "num_workers" in inf_cfg:
            cmd += [
                "--num_workers",
                str(inf_cfg["num_workers"]),
            ]

        if "cache" in inf_cfg:
            cmd += [
                "--cache",
                str(inf_cfg["cache"]),
            ]

        if inf_cfg.get("use_msa_server", False):
            cmd.append("--use_msa_server")

        if "msa_pairing_strategy" in inf_cfg:
            cmd += [
                "--msa_pairing_strategy",
                str(inf_cfg["msa_pairing_strategy"]),
            ]
        
        if use_template:
            cmd.append("--use_template")


        if inf_cfg.get("only_run_data_process", False):
            cmd.append("--only_run_data_process")

        runtime_params = {
            "run_id": run_id,
            "device": device,
            "model": inf_cfg.get("model", "v2-flash"),
            "precision": inf_cfg.get("precision"),
            "seed": inf_cfg.get("seed"),
            "num_workers": inf_cfg.get("num_workers"),
            "cache": inf_cfg.get("cache"),
        }

        sampling_params = {
            "num_diffusion_samples": num_diffusion_samples,
            "sampling_steps": sampling_steps,
            "recycling_iters": recycling_iters,
        }

        data_params = {
            "use_msa_server": inf_cfg.get("use_msa_server", False),
            "msa_pairing_strategy": inf_cfg.get("msa_pairing_strategy"),
            "use_template": use_template,
            "templates_available": len(available_templates),
            "templates_used": len(effective_templates),
            "only_run_data_process": inf_cfg.get("only_run_data_process", False),
        }

        runtime_params = {
            k: v for k, v in runtime_params.items()
            if v is not None
        }

        sampling_params = {
            k: v for k, v in sampling_params.items()
            if v is not None
        }

        data_params = {
            k: v for k, v in data_params.items()
            if v is not None
        }

        log_parameters(
            runtime_params,
            title=f"IntelliFold run={run_id} device={device} | Runtime"
        )

        log_parameters(
            sampling_params,
            title=f"IntelliFold run={run_id} device={device} | Sampling"
        )

        log_parameters(
            data_params,
            title=f"IntelliFold run={run_id} device={device} | Data/template options"
        )

        # runtime_params = {
        #     k: v for k, v in runtime_params.items()
        #     if v is not None
        # }

        log_parameters(
            runtime_params,
            title=f"IntelliFold run={run_id} device={device}"
        )

        subprocess.run(cmd, check=True, env=env)

        results = self._parse_results(out_dir)

        if not results:
            logger.warning(f"[WARNING] No IntelliFold samples returned for run {run_id}")

        return {
            "run_id": run_id,
            "device": device,
            "results": results,
            "tool": self.name,
        }

    def _parse_results(self, out_dir: Path):
        """
        Parse IntelliFold outputs.

        Handles two possible layouts:

        1. Boltz-like / IntelliFold observed layout:
        results/input/predictions/input/
            confidence_input_model_0.json
            input_model_0.cif
            confidence_input_model_1.json
            input_model_1.cif

        2. AF3-like layout:
        seed-42_sample-0/
            model.cif
            confidences.json
            summary_confidences.json

        The parser prefers confidence JSONs because they allow us to associate
        metrics with the correct structure.
        """

        samples = []

        # ============================================================
        # Case 1: IntelliFold observed layout
        #
        # results/input/predictions/input/
        #   input_seed-42_sample-0.cif
        #   input_seed-42_sample-0_confidences.json
        #   input_seed-42_sample-0_summary_confidences.json
        # ============================================================

        confidence_jsons = sorted(
            p for p in out_dir.rglob("*_confidences.json")
            if "_summary_confidences" not in p.name
        )

        if confidence_jsons:
            logger.info(
                f"[IntelliFold] Found {len(confidence_jsons)} confidence JSON files"
            )

            for confidence_json in confidence_jsons:
                parent_dir = confidence_json.parent

                # input_seed-42_sample-0_confidences.json
                # -> input_seed-42_sample-0
                model_stem = confidence_json.stem.replace(
                    "_confidences",
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
                        f"[IntelliFold] No structure file found for {confidence_json}"
                    )
                    continue

                summary_json = parent_dir / f"{model_stem}_summary_confidences.json"

                confidence = {}
                summary = {}

                if confidence_json.exists():
                    confidence = json.loads(confidence_json.read_text())

                if summary_json.exists():
                    summary = json.loads(summary_json.read_text())
                else:
                    logger.debug(
                        f"[IntelliFold] No summary confidence JSON found for {model_stem}"
                    )

                plddt = _first_present(
                    summary.get("plddt"),
                    summary.get("mean_plddt"),
                    summary.get("complex_plddt"),
                    confidence.get("plddt"),
                    confidence.get("mean_plddt"),
                    confidence.get("complex_plddt"),
                    _mean_if_list(confidence.get("atom_plddts")),
                    _mean_if_list(confidence.get("token_plddts")),
                )

                ptm = _first_present(
                    summary.get("ptm"),
                    summary.get("pTM"),
                    confidence.get("ptm"),
                    confidence.get("pTM"),
                )

                iptm = _first_present(
                    summary.get("iptm"),
                    summary.get("ipTM"),
                    confidence.get("iptm"),
                    confidence.get("ipTM"),
                )

                iplddt = _first_present(
                    summary.get("iplddt"),
                    summary.get("complex_iplddt"),
                    confidence.get("iplddt"),
                    confidence.get("complex_iplddt"),
                )

                ranking_score = _first_present(
                    summary.get("ranking_score"),
                    summary.get("ranking_confidence"),
                    summary.get("score"),
                    confidence.get("ranking_score"),
                    confidence.get("ranking_confidence"),
                    confidence.get("score"),
                )

                samples.append({
                    "structure": str(structure_file),
                    "plddt": plddt,
                    "ptm": ptm,
                    "iptm": iptm,
                    "iplddt": iplddt,
                    "ranking_score": ranking_score,
                })

                logger.debug(
                    f"[IntelliFold] Parsed {confidence_json.name} -> "
                    f"{structure_file.name} "
                    f"plddt={plddt} ptm={ptm} iptm={iptm} "
                    f"ranking_score={ranking_score}"
                )

            return samples

        # ============================================================
        # Case 2: AF3-like sample directories
        # ============================================================

        sample_dirs = sorted(
            [
                p for p in out_dir.rglob("*")
                if p.is_dir() and "sample" in p.name.lower()
            ]
        )

        if sample_dirs:
            logger.info(
                f"[IntelliFold] Found {len(sample_dirs)} sample directories"
            )

            for sample_dir in sample_dirs:
                structure_file = None

                for candidate in [
                    sample_dir / "model.cif",
                    *sorted(sample_dir.glob("*model*.cif")),
                    *sorted(sample_dir.glob("*.cif")),
                ]:
                    if candidate.exists():
                        structure_file = candidate
                        break

                if structure_file is None:
                    logger.warning(f"No CIF structure found in {sample_dir}")
                    continue

                summary_json = None
                confidence_json = None

                for candidate in [
                    sample_dir / "summary_confidences.json",
                    *sorted(sample_dir.glob("*summary*confidence*.json")),
                ]:
                    if candidate.exists():
                        summary_json = candidate
                        break

                for candidate in [
                    sample_dir / "confidences.json",
                    *sorted(sample_dir.glob("*confidence*.json")),
                ]:
                    if candidate.exists() and candidate != summary_json:
                        confidence_json = candidate
                        break

                summary = {}
                confidence = {}

                if summary_json and summary_json.exists():
                    summary = json.loads(summary_json.read_text())

                if confidence_json and confidence_json.exists():
                    confidence = json.loads(confidence_json.read_text())

                plddt = _first_present(
                    summary.get("plddt"),
                    summary.get("mean_plddt"),
                    summary.get("complex_plddt"),
                    _mean_if_list(confidence.get("atom_plddts")),
                    _mean_if_list(confidence.get("token_plddts")),
                )

                ptm = _first_present(
                    summary.get("ptm"),
                    summary.get("pTM"),
                )

                iptm = _first_present(
                    summary.get("iptm"),
                    summary.get("ipTM"),
                )

                iplddt = _first_present(
                    summary.get("iplddt"),
                    summary.get("complex_iplddt"),
                )

                ranking_score = _first_present(
                    summary.get("ranking_score"),
                    summary.get("ranking_confidence"),
                )

                samples.append({
                    "structure": str(structure_file),
                    "plddt": plddt,
                    "ptm": ptm,
                    "iptm": iptm,
                    "iplddt": iplddt,
                    "ranking_score": ranking_score,
                })

            return samples

        # ============================================================
        # Last-resort fallback: structures only, no metrics
        # ============================================================

        logger.warning(
            f"No IntelliFold confidence JSONs found under {out_dir}; "
            "falling back to recursive CIF discovery with empty metrics."
        )

        cif_files = sorted(out_dir.rglob("*.cif"))

        for cif_file in cif_files:
            samples.append({
                "structure": str(cif_file),
                "plddt": None,
                "ptm": None,
                "iptm": None,
                "iplddt": None,
                "ranking_score": None,
            })

        return samples