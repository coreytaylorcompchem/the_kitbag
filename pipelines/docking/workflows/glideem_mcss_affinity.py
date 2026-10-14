from __future__ import annotations

import csv
from pathlib import Path
import yaml

from backends import get_backend
from modules.glideem_tasks import analyse_correlations, extract_complexes
from pipeline.task_registry import get_task
from workflows import register_workflow


def read_manifest(path: Path):
    required = {"job_name", "receptor", "density_map", "reference_ligand"}
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"Empty GlideEM manifest: {path}")
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"Manifest missing columns: {sorted(missing)}")
    seen = set()
    for row in rows:
        if row["job_name"] in seen:
            raise ValueError(f"Duplicate job_name: {row['job_name']}")
        seen.add(row["job_name"])
        for field in required - {"job_name"}:
            p = Path(row[field]).expanduser()
            if not p.is_file():
                raise FileNotFoundError(f"{field} not found for {row['job_name']}: {p}")
    return rows


@register_workflow("glideem_mcss_affinity", description="GlideEM pose generation followed by MCSS docking and affinity correlation")
def run(config_path: str):
    with open(config_path) as handle:
        config = yaml.safe_load(handle)
    output = Path(config["output_dir"]).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    backend_cfg = dict(config["backend"])
    backend_name = backend_cfg.pop("name")
    backend = get_backend(backend_name, **backend_cfg)
    manifest_path = Path(config["inputs"]["glideem_manifest"]).expanduser().resolve()
    ligand_file = Path(config["inputs"]["activity_ligands_maegz"]).expanduser().resolve()
    if not ligand_file.is_file():
        raise FileNotFoundError(ligand_file)
    rows = read_manifest(manifest_path)

    prep_task = get_task("prepare_glideem_receptor")
    glideem_task = get_task("run_glideem")
    grid_task = get_task("build_glide_grid")
    mcss_task = get_task("dock_mcss")

    complex_jobs = []
    for row in rows:
        work = output / "glideem" / row["job_name"]
        work.mkdir(parents=True, exist_ok=True)
        prepared = prep_task(backend, None, config, job=row, work_dir=work)
        archive = glideem_task(backend, None, config, job=row, work_dir=work, prepared_receptor=prepared)
        complex_jobs.extend(extract_complexes(archive, output / "complexes", row["job_name"], config.get("force", False)))

    poseviewers = {}
    for job in complex_jobs:
        grid_task(backend, None, config, complex_job=job)
        poseviewers[job["name"]] = mcss_task(
            backend, None, config, complex_job=job, ligand_file=ligand_file, mcss_root=output / "mcss"
        )

    analyse_cfg = config.get("analysis", {})
    return analyse_correlations(
        complex_jobs, poseviewers, ligand_file,
        analyse_cfg["activity_column"], output / "analysis", analyse_cfg.get("min_pairs", 5)
    )
