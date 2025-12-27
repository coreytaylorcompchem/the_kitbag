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
        refinement: dict | None = None,
    ):
        """
        Parameters are injected by the StructureInferenceBackend scheduler.
        """

        # --- environment isolation ---
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = str(device)
        env["PYTHONHASHSEED"] = str(seed)

        output_dir.mkdir(parents=True, exist_ok=True)

        # --- input FASTA ---
        fasta = output_dir / "input.fasta"
        fasta.write_text(f">query\n{self.sequence}\n")

        # --- OpenFold output directory ---
        out_dir = output_dir / "results"
        out_dir.mkdir(exist_ok=True)

        # --- base OpenFold command ---
        cmd = [
            "openfold_run.py",
            "--fasta", str(fasta),
            "--output_dir", str(out_dir),
            "--seed", str(seed),
        ]

        # --- optional refinement / templating ---
        if refinement:
            reference = refinement.get("reference_structure")
            if reference:
                cmd += [
                    "--template_pdb", str(reference),
                    "--template_weight", str(refinement.get("weight", 0.5)),
                ]

            # region-based refinement is optional and tool-dependent
            if refinement.get("mode") == "region":
                for region in refinement.get("regions", []):
                    cmd += [
                        "--template_region",
                        f"{region['chain']}:{region['start']}-{region['end']}",
                    ]

        logger.info(
            f"[OpenFold] run={run_id} seed={seed} gpu={device}"
        )

        subprocess.run(cmd, check=True, env=env)

        # --- outputs ---
        structure = out_dir / "relaxed_model_1.pdb"
        metrics = self._parse_metrics(out_dir)

        return {
            "structure_path": str(structure),
            "metrics": metrics,
            "run_id": run_id,
            "seed": seed,
            "device": device,
        }

    def _parse_metrics(self, out_dir: Path):
        """
        Parse OpenFold confidence outputs.
        """
        import numpy as np

        plddt_path = out_dir / "plddt.npy"
        ptm_path = out_dir / "ptm.txt"
        iptm_path = out_dir / "iptm.txt"

        plddt = float(np.load(plddt_path).mean()) if plddt_path.exists() else None
        ptm = float(ptm_path.read_text()) if ptm_path.exists() else None
        iptm = float(iptm_path.read_text()) if iptm_path.exists() else None

        return {
            "plddt": plddt,
            "ptm": ptm,
            "iptm": iptm,
        }
