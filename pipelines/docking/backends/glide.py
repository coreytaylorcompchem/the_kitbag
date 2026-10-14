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
        path = self.schrodinger_root / name
        if not path.is_file():
            raise FileNotFoundError(f"Schrodinger executable not found: {path}")
        return path

    def run(self, cmd: list[Any], cwd: Path, log_path: Path) -> None:
        cwd.mkdir(parents=True, exist_ok=True)
        cmd = [str(x) for x in cmd]
        logger.info("Running: %s", shlex.join(cmd))
        with log_path.open("w") as log:
            proc = subprocess.run(cmd, cwd=str(cwd), stdout=log, stderr=subprocess.STDOUT, text=True)
        if proc.returncode:
            tail = "\n".join(log_path.read_text(errors="replace").splitlines()[-50:])
            raise RuntimeError(f"Command failed ({proc.returncode}): {shlex.join(cmd)}\n{tail}")

    def run_python(self, script: Path, args: list[Any], cwd: Path, log_path: Path) -> None:
        self.run([self.exe("run"), "python", script, *args], cwd, log_path)

    def run_glide(self, input_file: Path, cwd: Path, log_path: Path) -> None:
        self.run([self.exe("glide"), "-WAIT", "-HOST", self.host, input_file.name], cwd, log_path)

    def run_ligprep(self, input_file: Path, output_file: Path, cwd: Path, ph: float, pht: float) -> None:
        self.run([
            self.exe("ligprep"), "-ismi", input_file, "-omae", output_file,
            "-ph", ph, "-pht", pht, "-s", 1, "-WAIT"
        ], cwd, cwd / "ligprep.log")
        require_file(output_file, "LigPrep output")

    def run_glideem_template(self, values: Mapping[str, Any], cfg: Mapping[str, Any], cwd: Path) -> Path:
        """Run a site-specific GlideEM command template and validate its expected archive."""
        template = cfg.get("command_template")
        if not template:
            raise RuntimeError(
                "glideem.command_template is required. Copy the command emitted by your installed "
                "Maestro GlideEM panel and replace paths with {receptor}, {density_map}, {ligand}, "
                "{output}, {job_name}, {host}, and {schrodinger}."
            )
        rendered = template.format(**values)
        cmd = shlex.split(rendered)
        self.run(cmd, cwd, cwd / f"{values['job_name']}_glideem.log")
        output = Path(str(values["output"]))
        require_file(output, "GlideEM complex archive")
        return output


def require_file(path: Path, label: str) -> None:
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"{label} missing or empty: {path}")
