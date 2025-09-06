
import os
import subprocess
from backends.base import BaseBackend

class XTBBackend(BaseBackend):
    def optimise(self, xyz_file, config):
        charge = str(config.get("charge", 0))
        uhf = str(config.get("uhf", 0))
        gfn = str(config.get("gfn", 2))

        print(f"[xtb] Optimizing with GFN{gfn}, charge={charge}, uhf={uhf}")
        
        cmd = [
            "xtb", xyz_file,
            "--opt",
            "--gfn", gfn,
            "--chrg", charge,
            "--uhf", uhf
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print("[xtb] Error:", result.stderr)
            raise RuntimeError("xtb optimisation failed.")

        if not os.path.exists("xtbopt.xyz"):
            raise FileNotFoundError("xtbopt.xyz not found after optimisation.")

        with open("xtbopt.xyz") as f:
            optimised_xyz = f.read()

        return optimised_xyz