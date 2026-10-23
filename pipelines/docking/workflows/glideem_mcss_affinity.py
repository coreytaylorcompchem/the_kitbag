from __future__ import annotations

import csv
from pathlib import Path

import yaml

from concurrent.futures import ThreadPoolExecutor, as_completed

from backends import get_backend

from modules.utils.glideem_helpers import (
    extract_complexes,
)

from modules.utils.glideem_plotting import (
    analyse_correlations,
)

from pipeline.task_registry import get_task
from workflows import register_workflow

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def get_enabled_tasks(config: dict) -> set[str]:
    """
    Return the task names explicitly enabled for this workflow.

    The workflow section is required so workflow execution is
    declarative rather than implicitly hard-coded.
    """

    workflow = config.get(
        "workflow"
    )

    if workflow is None:
        raise ValueError(
            "Config must contain a 'workflow' section "
            "listing the tasks to run."
        )

    if not isinstance(
        workflow,
        list,
    ):
        raise ValueError(
            "'workflow' must be a list of task names."
        )

    enabled = {
        str(task).strip()
        for task in workflow
        if str(task).strip()
    }

    if not enabled:
        raise ValueError(
            "'workflow' contains no enabled tasks."
        )

    return enabled

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

    enabled_tasks = get_enabled_tasks(
        config
    )

    logger.info(
        "Enabled workflow tasks: %s",
        ", ".join(
            sorted(enabled_tasks)
        ),
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
    # Retrieve registered tasks
    # --------------------------------------------------

    ligprep_task = get_task(
        "prepare_activity_ligands"
    )

    extract_reference_task = get_task(
        "extract_reference_ligand"
    )

    prep_task = get_task(
        "prepare_glideem_receptor"
    )

    glideem_task = get_task(
        "run_glideem"
    )

    glideem_analysis_task = get_task(
        "analyse_glideem_poses"
    )

    glideem_summary_task = get_task(
        "summarise_glideem_ensembles"
    )

    glideem_campaign_task = get_task(
        "analyse_glideem_campaign_landscape"
    )

    grid_task = get_task(
        "build_glide_grid"
    )

    mcss_task = get_task(
        "dock_mcss"
    )

    mmgbsa_task = get_task(
        "calculate_mmgbsa"
    )

    pymol_task = get_task(
        "export_mcss_pymol"
    )

    available_tasks = {
        "prepare_activity_ligands": ligprep_task,
        "extract_reference_ligand": extract_reference_task,
        "prepare_glideem_receptor": prep_task,
        "run_glideem": glideem_task,
        "analyse_glideem_poses": glideem_analysis_task,
        "summarise_glideem_ensembles": glideem_summary_task,
        "analyse_glideem_campaign_landscape": (
            glideem_campaign_task
        ),
        "build_glide_grid": grid_task,
        "dock_mcss": mcss_task,
        "calculate_mmgbsa": mmgbsa_task,
        "export_mcss_pymol": pymol_task,
    }

    known_tasks = (
        set(available_tasks)
        | {"analyse_correlations"}
    )

    unknown_tasks = (
        enabled_tasks
        - known_tasks
    )

    if unknown_tasks:
        raise ValueError(
            "Unknown tasks in workflow: "
            + ", ".join(
                sorted(unknown_tasks)
            )
        )

    missing_tasks = [
        name
        for name in enabled_tasks
        if (
            name in available_tasks
            and available_tasks[name] is None
        )
    ]

    if missing_tasks:
        raise RuntimeError(
            "Enabled tasks are not registered: "
            + ", ".join(
                sorted(missing_tasks)
            )
        )

    # --------------------------------------------------
    # Determine downstream requirements
    # --------------------------------------------------

    pose_dependent_tasks = {
        "build_glide_grid",
        "dock_mcss",
        "calculate_mmgbsa",
        "analyse_correlations",
        "export_mcss_pymol",
    }

    needs_complex_jobs = bool(
        pose_dependent_tasks
        & enabled_tasks
    )

    needs_grids = bool(
        {
            "dock_mcss",
            "calculate_mmgbsa",
            "analyse_correlations",
            "export_mcss_pymol",
        }
        & enabled_tasks
    )

    needs_poseviewers = bool(
        {
            "calculate_mmgbsa",
            "analyse_correlations",
            "export_mcss_pymol",
        }
        & enabled_tasks
    )

    needs_activity_ligands = bool(
        {
            "prepare_activity_ligands",
            "dock_mcss",
            "analyse_correlations",
        }
        & enabled_tasks
    )

    # --------------------------------------------------
    # Activity ligand library
    # --------------------------------------------------

    ligand_file = (
        output
        / "ligprep"
        / "activity_ligands_ligprep.maegz"
    )

    if (
        "prepare_activity_ligands"
        in enabled_tasks
    ):
        logger.info(
            "Preparing experimental activity ligand library"
        )

        ligand_file = ligprep_task(
            backend,
            None,
            config,
        )

    ligand_file = (
        Path(ligand_file)
        .expanduser()
        .resolve()
    )

    if (
        needs_activity_ligands
        and not ligand_file.is_file()
    ):
        raise FileNotFoundError(
            "Prepared activity ligand library is required "
            "by the enabled workflow but was not found: "
            f"{ligand_file}. Enable "
            "'prepare_activity_ligands' or provide an "
            "existing cached output."
        )

    # --------------------------------------------------
    # Per-structure GlideEM stages
    # --------------------------------------------------

    complex_jobs = []
    failed_jobs = []
    glideem_analysis_results = []

    glideem_cfg = config.get(
        "glideem",
        {},
    )

    glideem_max_parallel = int(
        glideem_cfg.get(
            "max_parallel",
            1,
        )
    )

    glideem_max_parallel = max(
        1,
        min(
            glideem_max_parallel,
            len(rows),
        ),
    )

    logger.info(
        "Processing %d GlideEM structures "
        "with up to %d concurrent structure jobs",
        len(rows),
        glideem_max_parallel,
    )

    def run_glideem_structure(
        row,
    ):
        """
        Process one manifest structure through all enabled
        per-structure GlideEM stages.

        This function deliberately keeps dependent operations
        serial within a structure while allowing independent
        structures to run concurrently.
        """
        
        structure_failures = []

        job_name = row[
            "job_name"
        ]

        work = (
            output
            / "glideem"
            / job_name
        )

        work.mkdir(
            parents=True,
            exist_ok=True,
        )

        reference_ligand = (
            work
            / (
                f"{job_name}"
                "_reference_ligand.maegz"
            )
        )

        prepared = (
            work
            / f"{job_name}_prepared.maegz"
        )

        archive = (
            work
            / f"{job_name}_glideem-out.maegz"
        )

        glideem_csv = (
            work
            / f"{job_name}_glideem-out.csv"
        )

        analysis_result = None
        generated_complexes = []

        # ------------------------------------------
        # Extract reference ligand
        # ------------------------------------------

        if (
            "extract_reference_ligand"
            in enabled_tasks
        ):
            logger.info(
                "[%s] Extracting reference ligand",
                job_name,
            )

            reference_ligand = (
                extract_reference_task(
                    backend,
                    None,
                    config,
                    job=row,
                    work_dir=work,
                )
            )

            reference_ligand = Path(
                reference_ligand
            ).resolve()

            logger.info(
                "[%s] Reference ligand ready: %s",
                job_name,
                reference_ligand,
            )

        # ------------------------------------------
        # Prepare receptor
        # ------------------------------------------

        if (
            "prepare_glideem_receptor"
            in enabled_tasks
        ):
            logger.info(
                "[%s] Preparing receptor",
                job_name,
            )

            prepared = prep_task(
                backend,
                None,
                config,
                job=row,
                work_dir=work,
            )

            prepared = Path(
                prepared
            ).resolve()

            logger.info(
                "[%s] Prepared receptor ready: %s",
                job_name,
                prepared,
            )

        # ------------------------------------------
        # Run GlideEM
        # ------------------------------------------

        if (
            "run_glideem"
            in enabled_tasks
        ):
            if not prepared.is_file():
                raise FileNotFoundError(
                    f"[{job_name}] GlideEM requires "
                    f"a prepared receptor, but it was "
                    f"not found: {prepared}"
                )

            logger.info(
                "[%s] Starting GlideEM",
                job_name,
            )

            archive = glideem_task(
                backend,
                None,
                config,
                job=row,
                work_dir=work,
                prepared_receptor=prepared,
            )

            archive = Path(
                archive
            ).resolve()

            logger.info(
                "[%s] GlideEM completed: %s",
                job_name,
                archive,
            )

        # ------------------------------------------
        # Per-job GlideEM analysis
        # ------------------------------------------

        if (
            "analyse_glideem_poses"
            in enabled_tasks
        ):
            if not archive.is_file():
                raise FileNotFoundError(
                    f"[{job_name}] GlideEM pose analysis "
                    f"requires an existing GlideEM archive: "
                    f"{archive}"
                )

            if not glideem_csv.is_file():
                raise FileNotFoundError(
                    f"[{job_name}] GlideEM pose analysis "
                    f"requires the GlideEM CSV: "
                    f"{glideem_csv}"
                )

            if not reference_ligand.is_file():
                raise FileNotFoundError(
                    f"[{job_name}] GlideEM pose analysis "
                    f"requires the original reference "
                    f"ligand: {reference_ligand}"
                )

            logger.info(
                "[%s] Analysing GlideEM pose ensemble",
                job_name,
            )

        # ------------------------------------------
        # Load existing per-job GlideEM analysis
        # when campaign analysis requires it
        # ------------------------------------------

        elif (
            "analyse_glideem_campaign_landscape"
            in enabled_tasks
        ):
            analysis_output_dir = (
                output
                / "analysis"
                / "glideem"
                / job_name
            )

            metrics_file = (
                analysis_output_dir
                / "glideem_pose_metrics.csv"
            )

            clusters_file = (
                analysis_output_dir
                / "pose_clusters.csv"
            )

            pairwise_file = (
                analysis_output_dir
                / "pose_pairwise_rmsd.csv"
            )

            required_files = [
                archive,
                prepared,
                reference_ligand,
                metrics_file,
                clusters_file,
                pairwise_file,
            ]

            missing_files = [
                path
                for path in required_files
                if not Path(path).is_file()
            ]

            if missing_files:
                logger.warning(
                    "[%s] Existing GlideEM analysis is "
                    "incomplete; unavailable for campaign "
                    "analysis. Missing: %s",
                    job_name,
                    ", ".join(
                        str(path)
                        for path in missing_files
                    ),
                )

            else:
                reference_smiles = (
                    row.get(
                        "reference_smiles",
                        "",
                    )
                    or ""
                ).strip()

                analysis_result = {
                    "job_name": job_name,

                    "requested_nposes": int(
                        config.get(
                            "glideem",
                            {},
                        ).get(
                            "nposes",
                            5,
                        )
                    ),

                    "cluster_rmsd_cutoff": float(
                        config.get(
                            "glideem_analysis",
                            {},
                        ).get(
                            "cluster_rmsd_cutoff",
                            2.0,
                        )
                    ),

                    "metrics": metrics_file,
                    "clusters": clusters_file,
                    "pairwise_rmsd": pairwise_file,

                    "archive": archive,
                    "prepared_receptor": prepared,
                    "cryoem_reference": (
                        reference_ligand
                    ),

                    "reference_ligand_resname": (
                        row[
                            "reference_ligand"
                        ]
                        .strip()
                        .upper()
                    ),

                    "ligand_from_smiles": bool(
                        reference_smiles
                    ),
                }

                logger.info(
                    "[%s] Loaded existing GlideEM analysis "
                    "for campaign landscape analysis",
                    job_name,
                )


            try:
                analysis_result = (
                    glideem_analysis_task(
                        backend,
                        None,
                        config,
                        archive=archive,
                        glideem_csv=glideem_csv,
                        job=row,
                        cryoem_reference=(
                            reference_ligand
                        ),
                        prepared_receptor=prepared,
                        output_dir=(
                            output
                            / "analysis"
                            / "glideem"
                            / job_name
                        ),
                    )
                )

            except Exception as exc:
                logger.error(
                    "[%s] GlideEM pose analysis failed; "
                    "continuing with downstream processing: %s",
                    job_name,
                    exc,
                    exc_info=True,
                )

                structure_failures.append(
                    {
                        "job_name": job_name,
                        "stage": "glideem_analysis",
                        "error": str(exc),
                    }
                )

        # ------------------------------------------
        # Reconstruct downstream pose jobs
        # ------------------------------------------

        if needs_complex_jobs:
            if not archive.is_file():
                raise FileNotFoundError(
                    f"[{job_name}] Downstream workflow "
                    f"requires a GlideEM archive, but "
                    f"none exists: {archive}"
                )

            reference_smiles = (
                row.get(
                    "reference_smiles",
                    "",
                )
                or ""
            ).strip()

            logger.info(
                "[%s] Loading GlideEM pose jobs",
                job_name,
            )

            generated_complexes = (
                extract_complexes(
                    archive=archive,
                    out_root=(
                        output
                        / "complexes"
                    ),
                    source_name=job_name,
                    ligand_resname=(
                        row[
                            "reference_ligand"
                        ]
                    ),
                    force=config.get(
                        "force",
                        False,
                    ),
                    ligand_from_smiles=bool(
                        reference_smiles
                    ),
                )
            )

            if not generated_complexes:
                raise RuntimeError(
                    f"No usable GlideEM complexes "
                    f"were extracted from {archive}"
                )

            logger.info(
                "[%s] Loaded %d GlideEM pose jobs",
                job_name,
                len(
                    generated_complexes
                ),
            )

        return {
            "job_name": job_name,
            "complex_jobs": (
                generated_complexes
            ),
            "analysis_result": (
                analysis_result
            ),
            "failures": structure_failures,
        }

    # --------------------------------------------------
    # Execute structures in parallel
    # --------------------------------------------------

    with ThreadPoolExecutor(
        max_workers=glideem_max_parallel
    ) as executor:

        futures = {
            executor.submit(
                run_glideem_structure,
                row,
            ): row
            for row in rows
        }

        for future in as_completed(
            futures
        ):
            row = futures[
                future
            ]

            job_name = row[
                "job_name"
            ]

            try:
                result = future.result()

            except Exception as exc:
                logger.error(
                    "[%s] GlideEM structure processing "
                    "failed; this structure will be "
                    "skipped: %s",
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

            complex_jobs.extend(
                result[
                    "complex_jobs"
                ]
            )

            analysis_result = result[
                "analysis_result"
            ]

            if analysis_result is not None:
                glideem_analysis_results.append(
                    analysis_result
                )

            failed_jobs.extend(
                result[
                    "failures"
                ]
            )

    logger.info(
        "GlideEM structure processing finished: "
        "%d succeeded, %d failed",
        (
            len(rows)
            - sum(
                failure[
                    "stage"
                ] == "glideem"
                for failure in failed_jobs
            )
        ),
        sum(
            failure[
                "stage"
            ] == "glideem"
            for failure in failed_jobs
        ),
    )

    # --------------------------------------------------
    # Cross-structure GlideEM summary
    # --------------------------------------------------

    if (
        "summarise_glideem_ensembles"
        in enabled_tasks
    ):

        if glideem_analysis_results:

            logger.info(
                "Generating cross-structure "
                "GlideEM ensemble summary"
            )

            try:
                glideem_summary_task(
                    backend,
                    None,
                    config,
                    analysis_results=(
                        glideem_analysis_results
                    ),
                    output_dir=(
                        output
                        / "analysis"
                        / "glideem"
                    ),
                )

            except Exception as exc:
                logger.error(
                    "Cross-structure GlideEM ensemble "
                    "summary failed; continuing workflow: %s",
                    exc,
                    exc_info=True,
                )

                failed_jobs.append(
                    {
                        "job_name": (
                            "glideem_campaign_summary"
                        ),
                        "stage": (
                            "glideem_analysis"
                        ),
                        "error": str(exc),
                    }
                )

        else:
            logger.warning(
                "Cross-structure GlideEM summary was enabled, "
                "but no per-structure analysis results are "
                "available in this run; skipping."
            )

    # --------------------------------------------------
    # Campaign-level aligned 3D GlideEM landscape
    # --------------------------------------------------

    if (
        "analyse_glideem_campaign_landscape"
        in enabled_tasks
    ):
        if len(
            glideem_analysis_results
        ) >= 2:
            logger.info(
                "Generating aligned campaign-level "
                "GlideEM pose landscape"
            )

            try:
                glideem_campaign_task(
                    backend,
                    None,
                    config,
                    analysis_results=(
                        glideem_analysis_results
                    ),
                    output_dir=(
                        output
                        / "analysis"
                        / "glideem"
                        / "campaign_landscape"
                    ),
                )

            except Exception as exc:
                logger.error(
                    "Campaign GlideEM landscape analysis "
                    "failed; continuing workflow: %s",
                    exc,
                    exc_info=True,
                )

                failed_jobs.append(
                    {
                        "job_name": (
                            "glideem_campaign_landscape"
                        ),
                        "stage": (
                            "glideem_campaign_analysis"
                        ),
                        "error": str(
                            exc
                        ),
                    }
                )

        else:
            logger.warning(
                "Campaign GlideEM landscape analysis was "
                "enabled, but fewer than two successful "
                "per-structure analyses are available."
            )

    # --------------------------------------------------
    # If nothing downstream of GlideEM is enabled,
    # workflow can finish here.
    # --------------------------------------------------

    if not needs_complex_jobs:

        logger.info(
            "No pose-dependent downstream tasks enabled."
        )

        if failed_jobs:
            logger.warning(
                "Workflow completed with %d "
                "failed/skipped job(s).",
                len(failed_jobs),
            )

        return None

    if not complex_jobs:
        logger.error(
            "No usable GlideEM pose jobs are available "
            "for downstream processing."
        )
        return None

    # --------------------------------------------------
    # Glide grids
    # --------------------------------------------------

    grids = {}
    grid_jobs = []

    if (
        "build_glide_grid"
        in enabled_tasks
    ):

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
            "Building Glide grids for %d poses "
            "with up to %d concurrent jobs",
            len(complex_jobs),
            grid_max_parallel,
        )

        def run_grid_job(job):

            logger.info(
                f"[{job['name']}] "
                f"Building Glide grid"
            )

            grid = grid_task(
                backend,
                None,
                config,
                complex_job=job,
            )

            return job["name"], grid

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

            for future in as_completed(
                futures
            ):

                job = futures[
                    future
                ]

                try:
                    name, grid = (
                        future.result()
                    )

                except Exception as exc:
                    logger.error(
                        "[%s] Glide grid generation "
                        "failed; this pose will be "
                        "skipped: %s",
                        job["name"],
                        exc,
                        exc_info=True,
                    )

                    failed_jobs.append(
                        {
                            "job_name": (
                                job["name"]
                            ),
                            "stage": "grid",
                            "error": str(exc),
                        }
                    )

                    continue

                grids[name] = grid

                grid_jobs.append(
                    job
                )

                logger.debug(
                    f"[{name}] Glide grid ready: "
                    f"{grid}"
                )

    elif needs_grids:

        logger.info(
            "Grid generation is disabled; "
            "loading existing grids"
        )

        for job in complex_jobs:

            grid = Path(
                job["grid"]
            ).resolve()

            if grid.is_file():
                grids[
                    job["name"]
                ] = grid

                grid_jobs.append(
                    job
                )

            else:
                logger.warning(
                    "[%s] Existing grid not found; "
                    "pose will be unavailable downstream: %s",
                    job["name"],
                    grid,
                )

    # No downstream stage requires grids.
    else:
        grid_jobs = list(
            complex_jobs
        )

    if (
        needs_grids
        and not grid_jobs
    ):
        logger.error(
            "No usable Glide grids are available "
            "for downstream processing."
        )
        return None

    # --------------------------------------------------
    # MCSS docking / existing PoseViewers
    # --------------------------------------------------

    poseviewers = {}

    if "dock_mcss" in enabled_tasks:

        if not ligand_file.is_file():
            raise FileNotFoundError(
                "MCSS docking requires the prepared "
                f"activity ligand library: {ligand_file}"
            )

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
            "Starting MCSS docking against "
            "%d GlideEM poses with up to "
            "%d concurrent jobs",
            len(grid_jobs),
            max_parallel,
        )
        #debug
        mcss_host = mcss_cfg.get(
            "host",
            backend.host,
        )

        logger.debug(
            "MCSS concurrency: requested=%d, available_jobs=%d, "
            "effective_workers=%d, job_host=%s",
            int(
                mcss_cfg.get(
                    "max_parallel",
                    1,
                )
            ),
            len(grid_jobs),
            max_parallel,
            mcss_host,
        )

        def run_mcss_job(job):

            logger.info(
                f"[{job['name']}] "
                f"Starting MCSS docking"
            )

            pv = mcss_task(
                backend,
                None,
                config,
                complex_job=job,
                ligand_file=ligand_file,
                mcss_root=(
                    output
                    / "mcss"
                ),
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

            for future in as_completed(
                futures
            ):

                job = futures[
                    future
                ]

                try:
                    name, pv = (
                        future.result()
                    )

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
                            "job_name": (
                                job["name"]
                            ),
                            "stage": "mcss",
                            "error": str(exc),
                        }
                    )

                    continue

                poseviewers[
                    name
                ] = Path(pv)

                logger.info(
                    f"[{name}] MCSS docking "
                    f"completed: {pv}"
                )

    elif needs_poseviewers:

        logger.info(
            "MCSS docking is disabled; "
            "loading existing PoseViewers"
        )

        for job in grid_jobs:

            pv = (
                output
                / "mcss"
                / job["name"]
                / (
                    f"{job['name']}"
                    "_mcss_pv.maegz"
                )
            )

            if pv.is_file():
                poseviewers[
                    job["name"]
                ] = pv

            else:
                logger.warning(
                    "[%s] Existing MCSS PoseViewer "
                    "not found; unavailable for "
                    "downstream analysis: %s",
                    job["name"],
                    pv,
                )

    if (
        needs_poseviewers
        and not poseviewers
    ):
        logger.error(
            "No MCSS PoseViewers are available "
            "for downstream processing."
        )
        return None

    # --------------------------------------------------
    # Prime MM-GBSA rescoring
    # --------------------------------------------------

    mmgbsa_results = {}

    if (
        "calculate_mmgbsa"
        in enabled_tasks
    ):
        mmgbsa_cfg = config.get(
            "mmgbsa",
            {},
        )

        mmgbsa_max_parallel = int(
            mmgbsa_cfg.get(
                "max_parallel",
                1,
            )
        )

        mmgbsa_max_parallel = max(
            1,
            min(
                mmgbsa_max_parallel,
                len(poseviewers),
            ),
        )

        logger.info(
            "Starting Prime MM-GBSA for "
            "%d MCSS docking jobs with up to "
            "%d concurrent jobs",
            len(poseviewers),
            mmgbsa_max_parallel,
        )

        jobs_by_name = {
            job["name"]: job
            for job in grid_jobs
        }

        def run_mmgbsa_job(
            name,
            poseviewer,
        ):
            job = jobs_by_name[
                name
            ]

            logger.info(
                "[%s] Starting Prime MM-GBSA",
                name,
            )

            results_file = mmgbsa_task(
                backend,
                None,
                config,
                complex_job=job,
                poseviewer=poseviewer,
                mmgbsa_root=(
                    output
                    / "mmgbsa"
                ),
            )

            return (
                name,
                Path(
                    results_file
                ).resolve(),
            )

        with ThreadPoolExecutor(
            max_workers=mmgbsa_max_parallel
        ) as executor:
            futures = {
                executor.submit(
                    run_mmgbsa_job,
                    name,
                    poseviewer,
                ): name
                for name, poseviewer in (
                    poseviewers.items()
                )
            }

            for future in as_completed(
                futures
            ):
                name = futures[
                    future
                ]

                try:
                    (
                        result_name,
                        results_file,
                    ) = future.result()

                except Exception as exc:
                    logger.error(
                        "[%s] Prime MM-GBSA failed; "
                        "docking analysis can continue "
                        "without MM-GBSA values: %s",
                        name,
                        exc,
                        exc_info=True,
                    )

                    failed_jobs.append(
                        {
                            "job_name": name,
                            "stage": "mmgbsa",
                            "error": str(exc),
                        }
                    )

                    continue

                mmgbsa_results[
                    result_name
                ] = results_file

                logger.info(
                    "[%s] Prime MM-GBSA "
                    "completed: %s",
                    result_name,
                    results_file,
                )

    elif (
        "analyse_correlations"
        in enabled_tasks
    ):
        logger.info(
            "Prime MM-GBSA calculation is disabled; "
            "loading existing MM-GBSA results where present"
        )

        for job in grid_jobs:
            name = job[
                "name"
            ]

            results_file = (
                output
                / "mmgbsa"
                / name
                / (
                    f"{name}"
                    "_mmgbsa_results.csv"
                )
            )

            if results_file.is_file():
                mmgbsa_results[
                    name
                ] = results_file

    # --------------------------------------------------
    # Docking/activity correlation
    # --------------------------------------------------

    summary = None

    if (
        "analyse_correlations"
        in enabled_tasks
    ):

        analyse_cfg = config.get(
            "analysis",
            {}
        )

        if (
            "activity_column"
            not in analyse_cfg
        ):
            raise ValueError(
                "'analyse_correlations' is enabled, "
                "but analysis.activity_column is missing."
            )

        logger.info(
            "Starting docking/activity "
            "correlation analysis for %d poses",
            len(poseviewers),
        )

        try:
            summary = analyse_correlations(
                complex_jobs=complex_jobs,
                poseviewers=poseviewers,
                mmgbsa_results=mmgbsa_results,
                ligand_file=ligand_file,
                source_ligand_file=Path(
                    config[
                        "inputs"
                    ][
                        "activity_ligands"
                    ]
                ).expanduser().resolve(),
                activity_column=analyse_cfg[
                    "activity_column"
                ],
                output_dir=(
                    output
                    / "analysis"
                    / "mcss"
                ),
                min_pairs=analyse_cfg.get(
                    "min_pairs",
                    5,
                ),
            )

        except Exception as exc:
            logger.error(
                "Docking/activity correlation "
                "analysis failed: %s",
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

        else:
            logger.info(
                "Docking/activity correlation "
                "analysis completed"
            )

    # --------------------------------------------------
    # MCSS PyMOL export
    # --------------------------------------------------

    if (
        "export_mcss_pymol"
        in enabled_tasks
    ):

        pymol_cfg = config.get(
            "pymol",
            {}
        )

        top_n = int(
            pymol_cfg.get(
                "top_n",
                5,
            )
        )

        binding_site_cutoff = (
            pymol_cfg.get(
                "binding_site_cutoff"
            )
        )

        logger.info(
            "Exporting top %d MCSS poses "
            "per job for PyMOL visualisation",
            top_n,
        )

        jobs_by_name = {
            job["name"]: job
            for job in grid_jobs
        }

        for name, pv in (
            poseviewers.items()
        ):

            try:
                job = jobs_by_name[
                    name
                ]

                source_name = job[
                    "source_name"
                ]

                cryoem_reference = (
                    output
                    / "glideem"
                    / source_name
                    / (
                        f"{source_name}"
                        "_reference_ligand.maegz"
                    )
                )

                pymol_task(
                    backend,
                    None,
                    config,
                    complex_job=job,
                    poseviewer=pv,
                    cryoem_reference=(
                        cryoem_reference
                    ),
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
                        "stage": (
                            "pymol_export"
                        ),
                        "error": str(exc),
                    }
                )

    # --------------------------------------------------
    # Campaign completion
    # --------------------------------------------------

    if failed_jobs:

        logger.warning(
            "Workflow completed with %d "
            "failed/skipped job(s):",
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