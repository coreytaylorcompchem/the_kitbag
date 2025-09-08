import subprocess
from pathlib import Path
from typing import Tuple, Optional


class GninaBackend:
    def __init__(self, gnina_executable: str = "gnina", use_gpu: bool = True):
        """
        Initialize Gnina backend.

        Parameters:
        -----------
        gnina_executable : str
            Path to gnina binary.
        use_gpu : bool
            If False, adds --cpu flag to force CPU mode.
        """
        self.gnina_executable = gnina_executable
        self.use_gpu = use_gpu

    def dock(
        self,
        receptor_path: Path,
        ligand_path: Path,
        output_path: Path,
        center: Tuple[float, float, float],
        size: Tuple[float, float, float],
        extra_args: Optional[list] = None
    ):
        """
        Run docking using gnina.

        Parameters:
        -----------
        receptor_path : Path
            Path to the prepared receptor PDBQT file.
        ligand_path : Path
            Path to the prepared ligand PDBQT file.
        output_path : Path
            Path to write the docking result (usually SDF).
        center : Tuple[float, float, float]
            Docking box center coordinates (x, y, z).
        size : Tuple[float, float, float]
            Docking box dimensions (x, y, z).
        extra_args : list, optional
            Additional arguments to pass to gnina command.
        """

        cmd = [
            self.gnina_executable,
            "-r", str(receptor_path),
            "-l", str(ligand_path),
            "--center_x", str(center[0]),
            "--center_y", str(center[1]),
            "--center_z", str(center[2]),
            "--size_x", str(size[0]),
            "--size_y", str(size[1]),
            "--size_z", str(size[2]),
            "-o", str(output_path),
            "--num_modes", "20",
            "--exhaustiveness", "8"
        ]

        # Use --cpu only if GPU is explicitly disabled
        if not self.use_gpu:
            cmd.append("--cpu")

        if extra_args:
            cmd.extend(extra_args)

        print(f"[INFO] Running gnina command:\n{' '.join(cmd)}")

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"[ERROR] GNINA failed:\n{result.stderr}")
            raise RuntimeError(f"GNINA docking failed. Check logs above.")

        print(f"[INFO] Docking completed successfully. Output: {output_path}")
        return output_path
