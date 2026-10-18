from __future__ import annotations

import csv
from pathlib import Path

import yaml

from concurrent.futures import ThreadPoolExecutor, as_completed

from backends import get_backend
from modules.glideem_tasks import (
    analyse_correlations,
    extract_complexes,
)
from pipeline.task_registry import get_task
from workflows import register_workflow

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def read_manifest(path: Path):

    required = {
        "job_name",
        "receptor",
        "density_map",
        "reference_ligand",
    }

    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise ValueError(
            f"Empty GlideEM manifest: {path}"
        )

    missing = required - set(rows[0])

    if missing:
        raise ValueError(
            f"Manifest missing columns: {sorted(missing)}"
        )

    seen = set()

    for row in rows:

        if row["job_name"] in seen:
            raise ValueError(
                f"Duplicate job_name: {row['job_name']}"
            )

        seen.add(row["job_name"])

        for field in ("receptor", "density_map"):

            p = (
                Path(row[field])
                .expanduser()
                .resolve()
            )

            if not p.is_file():
                raise FileNotFoundError(
                    f"{field} not found for "
                    f"{row['job_name']}: {p}"
                )

        if not row["reference_ligand"].strip():
            raise ValueError(
                f"reference_ligand is empty for "
                f"{row['job_name']}"
            )
        
        # Optional chemistry override. Empty CSV cells arrive
        # as empty strings; normalize them here.
        row["reference_smiles"] = (
            row.get("reference_smiles", "")
            or ""
        ).strip()

    return rows


@register_workflow(
    "glideem_mcss_affinity",
    description=(
        "GlideEM pose generation followed by "
        "MCSS docking and affinity correlation"
    ),
)
def run(config_path: str):

    # --------------------------------------------------
    # Load configuration
    # --------------------------------------------------

    with open(config_path) as handle:
        config = yaml.safe_load(handle)

    output = (
        Path(config["output_dir"])
        .expanduser()
        .resolve()
    )

    output.mkdir(
        parents=True,
        exist_ok=True,
    )

    # --------------------------------------------------
    # Backend
    # --------------------------------------------------

    backend_cfg = dict(
        config["backend"]
    )

    backend_name = backend_cfg.pop(
        "name"
    )

    backend = get_backend(
        backend_name,
        **backend_cfg,
    )

    # --------------------------------------------------
    # Manifest
    # --------------------------------------------------

    manifest_path = (
        Path(
            config["inputs"]["glideem_manifest"]
        )
        .expanduser()
        .resolve()
    )

    if not manifest_path.is_file():
        raise FileNotFoundError(
            f"GlideEM manifest not found: "
            f"{manifest_path}"
        )

    rows = read_manifest(
        manifest_path
    )

    # --------------------------------------------------
    # Prepare experimental ligand library
    # --------------------------------------------------

    ligprep_task = get_task(
        "prepare_activity_ligands"
    )

    if ligprep_task is None:
        raise RuntimeError(
            "Task 'prepare_activity_ligands' "
            "is not registered."
        )

    ligand_file = ligprep_task(
        backend,
        None,
        config,
    )

    ligand_file = Path(
        ligand_file
    ).resolve()

    if not ligand_file.is_file():
        raise RuntimeError(
            f"LigPrep output was not created: "
            f"{ligand_file}"
        )

    # --------------------------------------------------
    # Retrieve remaining registered tasks
    # --------------------------------------------------

    extract_reference_task = get_task(
        "extract_reference_ligand"
    )
    
    prep_task = get_task(
        "prepare_glideem_receptor"
    )

    glideem_task = get_task(
        "run_glideem"
    )

    grid_task = get_task(
        "build_glide_grid"
    )

    mcss_task = get_task(
        "dock_mcss"
    )

    pymol_task = get_task(
        "export_mcss_pymol"
    )


    required_tasks = {
        "extract_reference_ligand": extract_reference_task,
        "prepare_glideem_receptor": prep_task,
        "run_glideem": glideem_task,
        "build_glide_grid": grid_task,
        "dock_mcss": mcss_task,
        "export_mcss_pymol": pymol_task,
    }

    missing_tasks = [
        name
        for name, task in required_tasks.items()
        if task is None
    ]

    if missing_tasks:
        raise RuntimeError(
            "Required tasks are not registered: "
            + ", ".join(missing_tasks)
        )

    # --------------------------------------------------
    # GlideEM
    # --------------------------------------------------

    complex_jobs = []
    failed_jobs = []

    for row in rows:

        job_name = row["job_name"]

        work = (
            output
            / "glideem"
            / job_name
        )

        work.mkdir(
            parents=True,
            exist_ok=True,
        )

        try:
            # ---------------------------------------------
            # Extract bound ligand from original PDB
            # ---------------------------------------------

            logger.info(
                f"[{job_name}] Extracting reference ligand"
            )

            reference_ligand = extract_reference_task(
                backend,
                None,
                config,
                job=row,
                work_dir=work,
            )

            logger.info(
                f"[{job_name}] Reference ligand ready: "
                f"{reference_ligand}"
            )

            # ---------------------------------------------
            # Prepare receptor
            # ---------------------------------------------

            logger.info(
                f"[{job_name}] Preparing receptor"
            )

            prepared = prep_task(
                backend,
                None,
                config,
                job=row,
                work_dir=work,
            )

            logger.info(
                f"[{job_name}] Prepared receptor ready: "
                f"{prepared}"
            )

            # ---------------------------------------------
            # GlideEM
            # ---------------------------------------------

            logger.info(
                f"[{job_name}] Starting GlideEM"
            )

            archive = glideem_task(
                backend,
                None,
                config,
                job=row,
                work_dir=work,
                prepared_receptor=prepared,
            )

            logger.info(
                f"[{job_name}] GlideEM completed: "
                f"{archive}"
            )

            # ---------------------------------------------
            # Extract GlideEM poses
            # ---------------------------------------------

            logger.info(
                f"[{job_name}] Extracting GlideEM poses"
            )

            reference_smiles = (
                row.get("reference_smiles", "")
                or ""
            ).strip()

            generated_complexes = extract_complexes(
                archive=archive,
                out_root=output / "complexes",
                source_name=job_name,
                ligand_resname=row["reference_ligand"],
                force=config.get("force", False),
                ligand_from_smiles=bool(reference_smiles),
            )

            if not generated_complexes:
                raise RuntimeError(
                    f"No usable GlideEM complexes were "
                    f"extracted from {archive}"
                )

            logger.info(
                f"[{job_name}] Extracted "
                f"{len(generated_complexes)} GlideEM poses"
            )

            complex_jobs.extend(
                generated_complexes
            )

        except Exception as exc:
            logger.error(
                "[%s] GlideEM processing failed; "
                "this structure will be skipped: %s",
                job_name,
                exc,
                exc_info=True,
            )

            failed_jobs.append(
                {
                    "job_name": job_name,
                    "stage": "glideem",
                    "error": str(exc),
                }
            )

            continue

    if not complex_jobs:
        logger.error(
            "No GlideEM jobs produced usable "
            "protein-ligand complexes. "
            "Nothing remains to process."
        )

        if failed_jobs:
            logger.error(
                "Failed/skipped jobs: %d",
                len(failed_jobs),
            )

        return None

    # --------------------------------------------------
    # Build Glide grids in parallel
    # --------------------------------------------------

    grid_cfg = config.get(
        "grid",
        {}
    )

    grid_max_parallel = int(
        grid_cfg.get(
            "max_parallel",
            1,
        )
    )

    grid_max_parallel = max(
        1,
        min(
            grid_max_parallel,
            len(complex_jobs),
        ),
    )

    logger.info(
        f"Building Glide grids for "
        f"{len(complex_jobs)} GlideEM poses "
        f"with up to {grid_max_parallel} concurrent jobs"
    )

    def run_grid_job(job):
        logger.info(
            f"[{job['name']}] Building Glide grid"
        )

        grid = grid_task(
            backend,
            None,
            config,
            complex_job=job,
        )

        return job["name"], grid

    grids = {}
    grid_jobs = []

    with ThreadPoolExecutor(
        max_workers=grid_max_parallel
    ) as executor:

        futures = {
            executor.submit(
                run_grid_job,
                job,
            ): job
            for job in complex_jobs
        }

        for future in as_completed(futures):

            job = futures[future]

            try:
                name, grid = future.result()

            except Exception as exc:
                logger.error(
                    "[%s] Glide grid generation failed; "
                    "this pose will be skipped: %s",
                    job["name"],
                    exc,
                    exc_info=True,
                )

                failed_jobs.append(
                    {
                        "job_name": job["name"],
                        "stage": "grid",
                        "error": str(exc),
                    }
                )

                continue

            grids[name] = grid

            grid_jobs.append(
                job
            )

            logger.info(
                f"[{name}] Glide grid ready: {grid}"
            )

    logger.info(
        "Glide grid generation finished: "
        "%d succeeded, %d failed",
        len(grid_jobs),
        len(complex_jobs) - len(grid_jobs),
    )

    if not grid_jobs:
        logger.error(
            "No Glide grids were generated successfully. "
            "MCSS docking cannot be performed."
        )
        return None


    # --------------------------------------------------
    # MCSS docking
    # --------------------------------------------------

    poseviewers = {}

    mcss_cfg = config.get(
        "mcss",
        {}
    )

    max_parallel = int(
        mcss_cfg.get(
            "max_parallel",
            1,
        )
    )

    max_parallel = max(
        1,
        min(
            max_parallel,
            len(grid_jobs),
        ),
    )

    logger.info(
        f"Starting MCSS docking against "
        f"{len(grid_jobs)} successfully gridded "
        f"GlideEM poses with up to "
        f"{max_parallel} concurrent jobs"
    )

    def run_mcss_job(job):
        logger.info(
            f"[{job['name']}] Starting MCSS docking"
        )

        pv = mcss_task(
            backend,
            None,
            config,
            complex_job=job,
            ligand_file=ligand_file,
            mcss_root=output / "mcss",
        )

        return job["name"], pv

    with ThreadPoolExecutor(
        max_workers=max_parallel
    ) as executor:

        futures = {
            executor.submit(
                run_mcss_job,
                job,
            ): job
            for job in grid_jobs
        }

        for future in as_completed(futures):

            job = futures[future]

            try:
                name, pv = future.result()

            except Exception as exc:
                logger.error(
                    "[%s] MCSS docking failed; "
                    "this pose will be skipped: %s",
                    job["name"],
                    exc,
                    exc_info=True,
                )

                failed_jobs.append(
                    {
                        "job_name": job["name"],
                        "stage": "mcss",
                        "error": str(exc),
                    }
                )

                continue

            poseviewers[name] = pv

            logger.info(
                f"[{name}] MCSS docking completed: "
                f"{pv}"
            )

    logger.info(
        "MCSS docking finished: "
        "%d succeeded, %d failed",
        len(poseviewers),
        len(grid_jobs) - len(poseviewers),
    )

    if not poseviewers:
        logger.error(
            "No MCSS docking jobs completed successfully. "
            "Skipping correlation analysis."
        )
        return None

    # --------------------------------------------------
    # Docking/activity correlation
    # --------------------------------------------------

    analyse_cfg = config.get(
        "analysis",
        {}
    )

    logger.info(
        "Starting docking/activity correlation analysis "
        f"for {len(poseviewers)} GlideEM poses"
    )

    try:
        summary = analyse_correlations(
            complex_jobs=complex_jobs,
            poseviewers=poseviewers,
            ligand_file=ligand_file,
            source_ligand_file=Path(
                config["inputs"]["activity_ligands"]
            ).expanduser().resolve(),
            activity_column=analyse_cfg["activity_column"],
            output_dir=output / "analysis",
            min_pairs=analyse_cfg.get("min_pairs", 5),
        )

    except Exception as exc:
        logger.error(
            "Docking/activity correlation analysis failed: %s",
            exc,
            exc_info=True,
        )

        failed_jobs.append(
            {
                "job_name": "analysis",
                "stage": "analysis",
                "error": str(exc),
            }
        )

        summary = None

    if summary is not None:
        logger.info(
            "Docking/activity correlation analysis completed"
        )

    # --------------------------------------------------
    # Output pymol script
    # --------------------------------------------------

    pymol_cfg = config.get(
        "pymol",
        {}
    )

    if pymol_cfg.get(
        "enabled",
        False,
    ):

        top_n = int(
            pymol_cfg.get(
                "top_n",
                5,
            )
        )

        logger.info(
            "Exporting top %d MCSS poses per job "
            "for PyMOL visualisation",
            top_n,
        )

        jobs_by_name = {
            job["name"]: job
            for job in grid_jobs
        }

        for name, pv in poseviewers.items():

            try:
                source_name = jobs_by_name[name]["source_name"]

                cryoem_reference = (
                    output
                    / "glideem"
                    / source_name
                    / f"{source_name}_reference_ligand.maegz"
                )

                pymol_task(
                    backend,
                    None,
                    config,
                    complex_job=jobs_by_name[name],
                    poseviewer=pv,
                    cryoem_reference=cryoem_reference,
                    top_n=top_n,
                )

            except Exception as exc:
                logger.error(
                    "[%s] PyMOL export failed; "
                    "continuing workflow: %s",
                    name,
                    exc,
                    exc_info=True,
                )

                failed_jobs.append(
                    {
                        "job_name": name,
                        "stage": "pymol_export",
                        "error": str(exc),
                    }
                )

    else:
        logger.info(
            "PyMOL export not enabled; skipping"
        )

        top_n = int(
            pymol_cfg.get(
                "top_n",
                5,
            )
        )

        jobs_by_name = {
            job["name"]: job
            for job in grid_jobs
        }

        for name, pv in poseviewers.items():

            try:
                source_name = jobs_by_name[name]["source_name"]

                cryoem_reference = (
                    output
                    / "glideem"
                    / source_name
                    / f"{source_name}_reference_ligand.maegz"
                )

                pymol_task(
                    backend,
                    None,
                    config,
                    complex_job=jobs_by_name[name],
                    poseviewer=pv,
                    cryoem_reference=cryoem_reference,
                    top_n=top_n,
                )

            except Exception as exc:
                logger.error(
                    "[%s] PyMOL export failed; "
                    "continuing workflow: %s",
                    name,
                    exc,
                    exc_info=True,
                )

                failed_jobs.append(
                    {
                        "job_name": name,
                        "stage": "pymol_export",
                        "error": str(exc),
                    }
                )

    if failed_jobs:
        logger.warning(
            "Workflow completed with %d failed/skipped job(s):",
            len(failed_jobs),
        )

        for failure in failed_jobs:
            logger.warning(
                "  [%s] stage=%s: %s",
                failure["job_name"],
                failure["stage"],
                failure["error"],
            )

    else:
        logger.info(
            "Workflow completed with no failed jobs."
        )

    return summary