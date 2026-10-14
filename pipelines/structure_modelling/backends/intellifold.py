import os
import yaml
import subprocess
import json
import statistics

from pathlib import Path

from backends.base import BaseStructureTool
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)


def validate_intellifold_msa_policy(
    sequences: dict,
    msas: dict,
    use_msa_server: bool,
):
    msas = msas or {}

    missing = [
        chain_id
        for chain_id in sequences
        if chain_id not in msas
    ]

    if missing and not use_msa_server:
        raise ValueError(
            "IntelliFold: use_msa_server=false, but no local MSA "
            "was provided for protein chains: "
            f"{', '.join(missing)}"
        )

    for chain_id, msa_path in msas.items():
        if chain_id not in sequences:
            logger.warning(
                f"[IntelliFold] Ignoring MSA for unknown chain "
                f"{chain_id}: {msa_path}"
            )
            continue

        msa_path = Path(msa_path).expanduser().resolve()

        if not msa_path.is_file():
            raise FileNotFoundError(
                f"IntelliFold MSA for chain {chain_id} "
                f"does not exist: {msa_path}"
            )

        if msa_path.stat().st_size == 0:
            raise ValueError(
                f"IntelliFold MSA for chain {chain_id} "
                f"is empty: {msa_path}"
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

def build_intellifold_command(
    yaml_file: Path,
    out_dir: Path,
    inf_cfg: dict,
):
    cmd = [
        "intellifold",
        "predict",
        str(yaml_file),
        "--out_dir",
        str(out_dir),
    ]

    # ------------------------------------------------------------
    # Generic pipeline names mapped to IntelliFold CLI names
    # ------------------------------------------------------------
    num_diffusion_samples = inf_cfg.get(
        "num_diffusion_samples",
        inf_cfg.get("diffusion_samples"),
    )

    recycling_iters = inf_cfg.get(
        "recycling_iters",
        inf_cfg.get("recycling_steps"),
    )

    sampling_steps = inf_cfg.get(
        "sampling_steps"
    )

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

    # ------------------------------------------------------------
    # Common scalar options
    # ------------------------------------------------------------
    scalar_options = {
        "seed": "--seed",
        "precision": "--precision",
        "num_workers": "--num_workers",
        "cache": "--cache",
        "msa_pairing_strategy": "--msa_pairing_strategy",
    }

    for config_key, cli_flag in scalar_options.items():
        value = inf_cfg.get(config_key)

        if value is not None:
            cmd += [
                cli_flag,
                str(value),
            ]

    # ------------------------------------------------------------
    # Boolean flags
    # ------------------------------------------------------------
    boolean_options = {
        "use_msa_server": "--use_msa_server",
        "use_template": "--use_template",
        "only_run_data_process": "--only_run_data_process",
        "no_pairing": "--no_pairing",
        "return_similar_seq": "--return_similar_seq",
        "override": "--override",
    }

    for config_key, cli_flag in boolean_options.items():
        if inf_cfg.get(config_key, False):
            cmd.append(cli_flag)

    # ------------------------------------------------------------
    # Optional escape hatch for version-specific IntelliFold flags
    # ------------------------------------------------------------
    extra_args = inf_cfg.get(
        "extra_args",
        []
    )

    if not isinstance(extra_args, list):
        raise TypeError(
            "structure_prediction.inference.extra_args "
            "must be a list"
        )

    cmd.extend(
        str(arg)
        for arg in extra_args
    )

    return cmd

def generate_intellifold_yaml(
    sequences: dict,
    ligands: list,
    templates: list,
    msas: dict,
    yaml_path: Path,
):
    """
    Generate IntelliFold input YAML.

    Protein entries contain:
      - chain ID
      - sequence
      - precomputed MSA path, when available

    Ligands are represented by SMILES.

    Templates retain the existing backend behaviour for now.
    """

    data = {
        "sequences": []
    }

    msas = msas or {}

    # ------------------------------------------------------------
    # Proteins
    # ------------------------------------------------------------
    for chain_id, seq in sequences.items():

        protein_entry = {
            "id": chain_id,
            "sequence": seq,
        }

        # Add precomputed local MSA if available
        if chain_id in msas:
            msa_path = Path(
                msas[chain_id]
            ).expanduser().resolve()

            if not msa_path.is_file():
                raise FileNotFoundError(
                    f"IntelliFold MSA for chain {chain_id} "
                    f"does not exist: {msa_path}"
                )

            if msa_path.stat().st_size == 0:
                raise ValueError(
                    f"IntelliFold MSA for chain {chain_id} "
                    f"is empty: {msa_path}"
                )

            protein_entry["msa"] = str(msa_path)

            logger.info(
                f"[IntelliFold] Chain {chain_id} "
                f"using local MSA: {msa_path}"
            )

        data["sequences"].append({
            "protein": protein_entry
        })

    # ------------------------------------------------------------
    # Ligands
    # ------------------------------------------------------------
    for i, smiles in enumerate(
        ligands or [],
        start=1,
    ):
        data["sequences"].append({
            "ligand": {
                "id": f"L{i}",
                "smiles": smiles,
            }
        })

    # ------------------------------------------------------------
    # Templates
    # ------------------------------------------------------------
    if templates:
        data["templates"] = []

        for template in templates:
            template_path = Path(
                template
            ).expanduser().resolve()

            if not template_path.is_file():
                raise FileNotFoundError(
                    f"IntelliFold template does not exist: "
                    f"{template_path}"
                )

            ext = template_path.suffix.lower()

            if ext in {".cif", ".mmcif"}:
                data["templates"].append({
                    "cif": str(template_path)
                })

            elif ext == ".pdb":
                data["templates"].append({
                    "pdb": str(template_path)
                })

            else:
                raise ValueError(
                    f"Unsupported IntelliFold template format: "
                    f"{template_path}"
                )

    # ------------------------------------------------------------
    # Final validation
    # ------------------------------------------------------------
    if not data["sequences"]:
        raise ValueError(
            "No valid sequences or ligands provided "
            "for IntelliFold input"
        )

    yaml_path.write_text(
        yaml.safe_dump(
            data,
            sort_keys=False,
        )
    )

    logger.info(
        f"[IntelliFold] Wrote input YAML: {yaml_path}"
    )

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
        msas: dict = None,
        row_constraints: list = None,
    ):
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = str(device)

        run_dir = output_dir / f"intellifold_run_{run_id}"
        run_dir.mkdir(parents=True, exist_ok=True)

        yaml_file = run_dir / "input.yaml"

        out_dir = run_dir / "results"

        sp_cfg = self.config.get(
            "structure_prediction",
            {}
        )

        inf_cfg = sp_cfg.get(
            "inference",
            {}
        )

        use_msa_server = bool(
            inf_cfg.get(
                "use_msa_server",
                False,
            )
        )

        use_template = bool(
            inf_cfg.get(
                "use_template",
                False,
            )
        )

        available_templates = templates or []

        if available_templates and not use_template:
            logger.info(
                f"[IntelliFold] Templates were provided for run {run_id}, "
                "but use_template=false; templates will not be included."
            )

            effective_templates = []

        else:
            effective_templates = available_templates

        validate_intellifold_msa_policy(
            sequences=sequences,
            msas=msas or {},
            use_msa_server=use_msa_server,
        )

        generate_intellifold_yaml(
            sequences=sequences,
            ligands=ligands or [],
            templates=effective_templates,
            msas=msas or {},
            yaml_path=yaml_file,
        )

        cmd = build_intellifold_command(
            yaml_file=yaml_file,
            out_dir=out_dir,
            inf_cfg=inf_cfg,
        )

        logger.info(
            "[IntelliFold] Command: "
            + " ".join(
                subprocess.list2cmdline([str(arg)])
                for arg in cmd
            )
        )

        msa_params = {
            "local_msa_chains": sorted(
                (msas or {}).keys()
            ),
            "local_msa_count": len(
                msas or {}
            ),
            "expected_protein_chains": len(
                sequences
            ),
            "use_msa_server": use_msa_server,
            "msa_pairing_strategy": inf_cfg.get(
                "msa_pairing_strategy"
            ),
        }

        log_parameters(
            msa_params,
            title=(
                f"IntelliFold run={run_id} device={device} "
                "| MSA"
            ),
        )

        subprocess.run(
            cmd,
            check=True,
            env=env,
        )

        results = self._parse_results(out_dir)

        if not results:
            logger.warning(
                f"[WARNING] No IntelliFold samples returned "
                f"for run {run_id}"
            )

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