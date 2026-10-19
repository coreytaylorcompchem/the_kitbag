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
                stdin=subprocess.DEVNULL,
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

    def run_glide(
        self,
        input_file: Path,
        cwd: Path,
        log_path: Path,
        host: str | None = None,
    ) -> None:

        effective_host = (
            host or self.host
        )

        self.run(
            [
                self.exe("glide"),
                "-WAIT",
                "-HOST",
                effective_host,
                input_file.name,
            ],
            cwd,
            log_path,
        )

    def run_glideem(
        self,
        receptor: Path,
        density_map: Path,
        job_name: str,
        cwd: Path,
        mode: str = "exhaustive",
        nposes: int = 5,
        ligand_asl: str | None = None,
        ligand_smiles: str | None = None,
        docking_site_xyz: tuple[
            float,
            float,
            float,
        ] | None = None,
        host: str | None = None,
    ) -> None:
        """
        Run GlideEM using either:

        1. a prepared ligand already present in the receptor complex,
        selected by ligand_asl; or

        2. an authoritative ligand SMILES plus an experimental
        docking-site center.
        """

        if bool(ligand_asl) == bool(ligand_smiles):
            raise ValueError(
                "Exactly one of ligand_asl or ligand_smiles "
                "must be supplied."
            )

        cmd = [
            self.exe("run"),
            "-FROM",
            "psp",
            "glideEM.py",
            receptor,
            "-map",
            density_map,
            "-mode",
            mode,
            "-nposes",
            nposes,
        ]

        if ligand_smiles:

            if docking_site_xyz is None:
                raise ValueError(
                    "docking_site_xyz is required when "
                    "ligand_smiles is supplied."
                )

            x, y, z = docking_site_xyz

            cmd.extend(
                [
                    "-ligand_smiles",
                    ligand_smiles,
                    "-docking_site_xyz",
                    x,
                    y,
                    z,
                ]
            )

        else:
            cmd.extend(
                [
                    "-ligand_asl",
                    ligand_asl,
                ]
            )

        cmd.extend(
            [
                "-jobname",
                job_name,
                "-HOST",
                self.host,
                "-WAIT",
            ]
        )

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
