import os
import subprocess
import logging
from pathlib import Path

from backends.base import BaseStructureTool

logger = logging.getLogger(__name__)


class BoltzBackend(BaseStructureTool):
    name = "boltz"

    def run(
        self,
        run_id: int,
        seed: int,
        device: int,
        output_dir: Path,
        refinement: dict | None = None,
    ):
        """
        Run a single Boltz inference job.
        """

        # --- environment isolation ---
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = str(device)
        env["PYTHONHASHSEED"] = str(seed)

        output_dir.mkdir(parents=True, exist_ok=True)

        # --- input FASTA ---
        fasta = output_dir / "input.fasta"
        fasta.write_text(f">query\n{self.sequence}\n")

        # --- output directory ---
        out_dir = output_dir / "results"
        out_dir.mkdir(exist_ok=True)

        # --- base Boltz command ---
        cmd = [
            "boltz",
            "predict",
            str(fasta),
            "--outdir", str(out_dir),
            "--seed", str(seed),
        ]

        # --- optional refinement ---
        if refinement:
            reference = refinement.get("reference_structure")
            if reference:
                cmd += [
                    "--template", str(reference),
                    "--template_weight", str(refinement.get("weight", 0.5)),
                ]

            if refinement.get("mode") == "region":
                for region in refinement.get("regions", []):
                    cmd += [
                        "--template_region",
                        f"{region['chain']}:{region['start']}-{region['end']}",
                    ]

        logger.info(
            f"[Boltz] run={run_id} seed={seed} gpu={device}"
        )

        subprocess.run(cmd, check=True, env=env)

        structure = out_dir / "ranked_0.pdb"
        metrics = self._parse_metrics(out_dir)

        return {
            "structure_path": str(structure),
            "metrics": metrics,
            "run_id": run_id,
            "seed": seed,
            "device": device,
        }

    def _parse_metrics(self, out_dir: Path):
        import json

        confidence_file = out_dir / "confidence.json"
        if not confidence_file.exists():
            logger.warning("[Boltz] confidence.json not found")
            return {"plddt": None, "ptm": None, "iptm": None}

        metrics = json.loads(confidence_file.read_text())

        return {
            "plddt": metrics.get("mean_plddt"),
            "ptm": metrics.get("ptm"),
            "iptm": metrics.get("iptm"),
        }
