import os
import subprocess
import shutil
import tempfile
from pathlib import Path
from backends.base import BaseBackend

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

class XTBBackend(BaseBackend):
    def optimise(self, xyz_file, config):
        charge = str(config.get("charge", 0))
        uhf = str(config.get("uhf", 0))
        gfn = str(config.get("gfn", 2))

        logger.info(f"Optimizing {xyz_file} with GFN{gfn}, uhf={uhf}")

        # Create a temporary directory for the XTB run
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            input_file = tmp_path / "input.xyz"
            shutil.copy(xyz_file, input_file)

            cmd = [
                "xtb", str(input_file),
                "--opt",
                "--gfn", gfn,
                # "--chrg", charge,
                "--uhf", uhf,
                "-alpb", "water"
            ]

            result = subprocess.run(cmd, cwd=tmp_path, capture_output=True, text=True)

            if result.returncode != 0:
                logger.error(f"xtb error: {result.stderr}")
                raise RuntimeError("xtb optimisation failed.")

            opt_file = tmp_path / "xtbopt.xyz"
            out_file = tmp_path / "xtb.out"

            if not opt_file.exists():
                raise FileNotFoundError("xtbopt.xyz not found after optimisation.")

            with open(opt_file) as f:
                optimised_xyz = f.read()

        
            energy = None
            if out_file.exists():
                with open(out_file) as f:
                    for line in f:
                        if "TOTAL ENERGY" in line:
                            energy = float(line.strip().split()[-1])
                            break

            return optimised_xyz, energy
