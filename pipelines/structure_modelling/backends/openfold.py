import os
import subprocess
from pathlib import Path
import logging

from backends.base import BaseStructureTool

logger = logging.getLogger(__name__)

class OpenFoldBackend(BaseStructureTool):
    name = "openfold"

    def run(
        self,
        run_id: int,
        seed: int,
        device: int,
        output_dir: Path,
        refinement: dict | None = None
    ):
        """
        Run OpenFold inference on the given sequence.

        Parameters:
        - run_id: unique run identifier (used for output dir)
        - seed: random seed for reproducibility
        - device: CUDA device ID
        - output_dir: base output directory
        - refinement: optional dict with {"pdb_path": ...} for local refinement
        """

        # --- set environment ---
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = str(device)
        env["PYTHONHASHSEED"] = str(seed)

        # --- prepare input fasta (supports multi-chain sequences) ---
        fasta = self.prepare_input(output_dir)

        # --- create per-run output dir ---
        out_dir = output_dir / f"run_{run_id}"
        out_dir.mkdir(exist_ok=True)

        # --- construct command ---
        cmd = [
            "openfold_run.py",
            "--fasta", str(fasta),
            "--output_dir", str(out_dir),
        ]

        if refinement:
            cmd.extend(["--refine_pdb", refinement["pdb_path"]])

        # --- execute ---
        subprocess.run(cmd, check=True, env=env)

        # --- parse results ---
        structure = out_dir / "relaxed_model_1.pdb"
        metrics = self._parse_metrics(out_dir)

        return {
            "run_id": run_id,
            "structure_path": str(structure),
            "metrics": metrics,
            "tool": self.name,
        }

    def _parse_metrics(self, out_dir: Path):
        """Parse PLDDT and pTM metrics safely."""
        import numpy as np

        plddt_file = out_dir / "plddt.npy"
        ptm_file = out_dir / "ptm.txt"

        plddt = float(np.load(plddt_file).mean()) if plddt_file.exists() else None
        ptm = float(ptm_file.read_text().strip()) if ptm_file.exists() else None

        return {
            "plddt": plddt,
            "ptm": ptm,
            "iptm": None,
        }

    def run_multi(
        self,
        seeds: list[int],
        device: int,
        output_dir: Path,
        refinement: dict | None = None
    ):
        """
        Run multiple OpenFold inferences with different seeds.
        Stores all results in self.cache['ranking']['results'].
        """
        all_results = []
        for run_id, seed in enumerate(seeds):
            result = self.run(run_id, seed, device, output_dir, refinement)
            all_results.append(result)

        # cache results for consensus / ranking
        self.cache = {"ranking": {"results": all_results}}
        return all_results