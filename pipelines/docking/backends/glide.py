from __future__ import annotations

import os
import shlex
import subprocess
from pathlib import Path
from typing import Any, Mapping

from backends.base import BaseBackend
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

class SchrodingerGlideEMBackend(BaseBackend):
    """Schrodinger backend for GlideEM, grid generation, LigPrep and Glide MCSS."""

    supported_tasks = [
        "prepare_glideem_receptor",
        "run_glideem",
        "extract_glideem_complexes",
        "prepare_activity_ligands",
        "build_glide_grid",
        "dock_mcss",
    ]

    def __init__(self, schrodinger_root: str | None = None, host: str = "localhost:1", **kwargs):
        root = schrodinger_root or os.environ.get("SCHRODINGER")
        if not root:
            raise RuntimeError("Set backend.schrodinger_root or the SCHRODINGER environment variable")
        self.schrodinger_root = Path(root).expanduser().resolve()
        self.host = host
        super().__init__(executable_path=str(self.schrodinger_root / "glide"), use_gpu=False)

    def exe(self, name: str) -> Path:
        """
        Resolve a Schrodinger executable.

        Schrodinger exposes some applications directly under
        $SCHRODINGER and others under $SCHRODINGER/utilities.
        """

        candidates = [
            self.schrodinger_root / name,
            self.schrodinger_root / "utilities" / name,
        ]

        for path in candidates:
            if path.is_file():
                return path

        searched = "\n".join(
            f"  - {path}"
            for path in candidates
        )

        raise FileNotFoundError(
            f"Schrodinger executable '{name}' not found.\n"
            f"Searched:\n{searched}"
        )

    def run(
        self,
        cmd: list[Any],
        cwd: Path,
        log_path: Path,
    ) -> None:

        cwd.mkdir(
            parents=True,
            exist_ok=True
        )

        cmd = [
            str(x)
            for x in cmd
        ]

        command_text = shlex.join(cmd)

        logger.info(
            "Running: %s",
            command_text,
        )

        logger.info(
            "Log: %s",
            log_path,
        )

        with log_path.open("w") as log:

            proc = subprocess.Popen(
                cmd,
                cwd=str(cwd),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )

            assert proc.stdout is not None

            for line in proc.stdout:

                line = line.rstrip()

                # Always preserve raw subprocess output
                # in the job log.
                log.write(line + "\n")
                log.flush()

                # Also expose it while the workflow runs.
                if line:
                    logger.info(
                        "[Schrodinger] %s",
                        line,
                    )

            returncode = proc.wait()

        if returncode:

            tail = "\n".join(
                log_path
                .read_text(errors="replace")
                .splitlines()[-50:]
            )

            raise RuntimeError(
                f"Command failed ({returncode}): "
                f"{command_text}\n"
                f"{tail}"
            )

        logger.info(
            "Command completed successfully."
        )

    def run_python(self, script: Path, args: list[Any], cwd: Path, log_path: Path) -> None:
        self.run([self.exe("run"), "python3", script, *args], cwd, log_path)

    def run_glide(self, input_file: Path, cwd: Path, log_path: Path) -> None:
        self.run([self.exe("glide"), "-WAIT", "-HOST", self.host, input_file.name], cwd, log_path)

    # def run_ligprep(self, input_file: Path, output_file: Path, cwd: Path, ph: float, pht: float) -> None:
    #     self.run([
    #         self.exe("ligprep"), "-ismi", input_file, "-omae", output_file,
    #         "-ph", ph, "-pht", pht, "-s", 1, "-WAIT"
    #     ], cwd, cwd / "ligprep.log")
    #     require_file(output_file, "LigPrep output")

    def run_glideem(
        self,
        receptor: Path,
        density_map: Path,
        ligand_asl: str,
        job_name: str,
        cwd: Path,
        mode: str = "exhaustive",
        nposes: int = 5,
    ) -> None:
        """
        Run GlideEM using a ligand already present in the prepared
        receptor complex.

        The ligand is selected from the prepared complex with an ASL,
        for example:

            res.ptype LIG
        """

        cmd = [
            self.exe("run"),
            "-FROM",
            "psp",
            "glideEM.py",

            receptor,

            "-map",
            density_map,

            "-ligand_asl",
            ligand_asl,

            "-mode",
            mode,

            "-nposes",
            nposes,

            "-jobname",
            job_name,

            "-HOST",
            self.host,

            "-WAIT",
        ]

        # This is Kitbag's launcher/captured-console log.
        # Schrodinger Job Control separately creates:
        #
        #     <job_name>.log
        #
        launcher_log = (
            cwd /
            f"{job_name}_launcher.log"
        )

        self.run(
            cmd,
            cwd,
            launcher_log,
        )

    def dock(self, ligand: dict, config: dict, **kwargs):
        raise NotImplementedError(
            "SchrodingerGlideEMBackend does not support generic dock(). "
            "Use the dock_mcss task."
        )

def require_file(path: Path, label: str) -> None:
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"{label} missing or empty: {path}")
