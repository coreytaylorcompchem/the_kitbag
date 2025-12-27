import os
import subprocess
import logging
from pathlib import Path
from backends.base import BaseStructureTool

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=True, simple_format=True)

class ChaiBackend(BaseStructureTool):
    name = "chai"

    def run(
        self,
        run_id: int,
        seed: int,
        device: int,
        output_dir: Path,
        refinement: dict | None = None
    ):
        """Run a single Chai inference job with optional refinement."""
        # --- environment isolation ---
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = str(device)
        env["PYTHONHASHSEED"] = str(seed)

        # --- per-run output dir ---
        run_dir = output_dir / f"run_{run_id}"
        run_dir.mkdir(parents=True, exist_ok=True)

        # --- input FASTA ---
        fasta = run_dir / "input.fasta"
        fasta.write_text(f">query\n{self.sequence}\n")

        # --- results dir ---
        out_dir = run_dir / "results"
        out_dir.mkdir(exist_ok=True)

        # --- base Chai command ---
        cmd = [
            "chai_predict",
            "--fasta", str(fasta),
            "--out", str(out_dir),
            "--seed", str(seed),
        ]

        # --- optional refinement / conditioning ---
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

        logger.info(f"[Chai] run={run_id} seed={seed} gpu={device}")
        subprocess.run(cmd, check=True, env=env)

        structure = out_dir / "model.pdb"
        metrics = self._parse_metrics(out_dir)

        return {
            "run_id": run_id,
            "seed": seed,
            "device": device,
            "structure_path": str(structure),
            "metrics": metrics,
            "tool": self.name,
        }

    def run_multi(
        self,
        seeds: list[int],
        device: int,
        output_dir: Path,
        refinement: dict | None = None
    ):
        all_results = []
        for run_id, seed in enumerate(seeds):
            result = self.run(run_id, seed, device, output_dir, refinement)
            all_results.append(result)

        self.cache = {"ranking": {"results": all_results}}
        return all_results

    def _parse_metrics(self, out_dir: Path):
        import json

        metrics_file = out_dir / "metrics.json"
        if not metrics_file.exists():
            logger.warning("[Chai] metrics.json not found")
            return {"plddt": None, "ptm": None, "iptm": None}

        data = json.loads(metrics_file.read_text())
        return {
            "plddt": data.get("plddt"),
            "ptm": data.get("ptm"),
            "iptm": data.get("iptm"),
        }
