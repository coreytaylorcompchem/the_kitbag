from __future__ import annotations

import re
import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.colors import (
            ListedColormap,
            BoundaryNorm,
        )
from matplotlib.patches import Patch

from schrodinger import structure

from scipy.cluster.hierarchy import (
    dendrogram,
    fcluster,
    linkage,
)
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

from modules.utils.glideem_helpers import(
    prepare_activity_csv_for_ligprep,
    annotate_ligprep_csv_metadata,
    heavy_atoms,
    glide_score,
    compound_id,
    annotate_ligprep_compound_ids,
    find_glideem_ligand_indices,
    ligand_heavy_atom_coordinates,
    get_reference_ligand_atoms,
    get_reference_ligand_center,
    coordinate_frame_rmsd,
    export_glideem_pymol_ensemble,
)

from modules.utils.glideem_plotting import (
    save_cluster_composition,
    save_cluster_dendrogram,
    save_clustered_score_vs_rmsd,
    save_density_rmsd_small_multiples,
    save_pose_ambiguity_map,
    save_reference_rmsd_distributions,
    save_rmsd_heatmap,
)

from modules.utils.pymol_export import (
    export_mcss_pose_bundle,
)

from pipeline.logger import setup_logger
from pipeline.task_registry import register_task

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_task(
    "prepare_activity_ligands",
    category="Ligand preparation",
    description="Prepare experimental activity ligands using Schrodinger LigPrep."
)
def prepare_activity_ligands(
    backend,
    ligand,
    config,
    **kwargs
):
    """
    Prepare the experimental ligand library for Glide docking.

    Input:
        config["inputs"]["activity_ligands"]

    Output:
        <output_dir>/ligprep/activity_ligands_ligprep.maegz
    """

    input_path = (
        Path(config["inputs"]["activity_ligands"])
        .expanduser()
        .resolve()
    )

    if not input_path.is_file():
        raise FileNotFoundError(
            f"Activity ligand input not found: {input_path}"
        )

    input_suffix = (
        input_path.suffix
        .lower()
    )

    if input_suffix not in {
        ".sdf",
        ".sd",
        ".csv",
    }:
        raise ValueError(
            "activity_ligands must be an SDF or CSV file; "
            f"received: {input_path}"
        )

    output_root = Path(config["output_dir"]).expanduser().resolve()

    output_dir = output_root / "ligprep"
    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    ligprep_cfg = config.get(
        "ligprep",
        {}
    )

    activity_column = (
        config
        .get("analysis", {})
        .get("activity_column")
    )

    ligprep_input = input_path
    csv_metadata = None

    if input_suffix == ".csv":

        if not activity_column:
            raise ValueError(
                "CSV activity ligand input requires "
                "analysis.activity_column in the YAML."
            )

        smiles_column = ligprep_cfg.get(
            "smiles_column",
            "SMILES",
        )

        id_column = ligprep_cfg.get(
            "id_column",
            "ID",
        )

        (
            ligprep_input,
            csv_metadata,
        ) = prepare_activity_csv_for_ligprep(
            input_csv=input_path,
            work_dir=output_dir,
            smiles_column=smiles_column,
            id_column=id_column,
            activity_column=activity_column,
        )

        logger.info(
            "Converted activity CSV to LigPrep "
            "SMILES input: %s",
            ligprep_input,
        )

    output_file = (
        output_dir /
        "activity_ligands_ligprep.maegz"
    )

    # Reuse existing output unless explicitly forcing rerun

    source_record = (
        output_dir
        / "activity_ligands_source.txt"
    )

    can_reuse = (
        output_file.is_file()
        and output_file.stat().st_size > 0
        and source_record.is_file()
        and source_record.read_text().strip()
        == str(input_path)
        and not config.get(
            "force",
            False,
        )
    )

    if can_reuse:
        logger.info(
            "Using existing LigPrep output: %s",
            output_file,
        )

        return output_file

        
    # if (
    #     output_file.is_file()
    #     and output_file.stat().st_size > 0
    #     and not config.get("force", False)
    # ):
    #     logger.info(
    #         f"Using existing LigPrep output: {output_file}"
    #     )

    #     return output_file

    ligprep_cfg = config.get(
        "ligprep",
        {}
    )

    pH = ligprep_cfg.get(
        "pH",
        7.4
    )

    pH_tolerance = ligprep_cfg.get(
        "pH_tolerance",
        2.0
    )

    max_states = ligprep_cfg.get(
        "max_states",
        1
    )

    if input_suffix == ".csv":
        input_args = [
            "-ismi",
            str(ligprep_input),
        ]

    else:
        input_args = [
            "-isd",
            str(ligprep_input),
        ]

    cmd = [
        backend.exe("ligprep"),
        *input_args,
        "-omae",
        str(output_file),
        "-ph",
        str(pH),
        "-pht",
        str(pH_tolerance),
        "-s",
        str(max_states),
        "-WAIT",
    ]

    log_file = (
        output_dir /
        "activity_ligands_ligprep.log"
    )

    logger.info(
        f"Preparing activity ligand library with LigPrep: "
        f"{input_path.name}"
    )

    backend.run(
        cmd,
        output_dir,
        log_file
    )

    if (
        not output_file.is_file()
        or output_file.stat().st_size == 0
    ):
        raise RuntimeError(
            f"LigPrep did not produce expected output: "
            f"{output_file}"
        )

    if input_suffix == ".csv":

        annotate_ligprep_csv_metadata(
            ligprep_file=output_file,
            metadata_by_id=csv_metadata,
            activity_column=activity_column,
        )

    else:

        annotate_ligprep_compound_ids(
            input_sdf=input_path,
            ligprep_file=output_file,
        )

    logger.info(
        f"LigPrep completed: {output_file}"
    )

    return output_file

@register_task(
    "extract_reference_ligand",
    category="GlideEM",
    description=(
        "Extract the specified reference ligand "
        "from the input protein-ligand PDB."
    ),
)

def extract_reference_ligand(
    backend,
    ligand,
    config,
    **kwargs,
):
    """
    Extract a bound ligand from an input PDB using
    the residue name supplied in the GlideEM manifest.

    Example manifest value:

        reference_ligand = LIG

    Output:

        <work_dir>/<job_name>_reference_ligand.maegz
    """

    job = kwargs["job"]
    work_dir = Path(kwargs["work_dir"])

    pdb_path = (
        Path(job["receptor"])
        .expanduser()
        .resolve()
    )

    ligand_resname = (
        job["reference_ligand"]
        .strip()
        .upper()
    )

    if not pdb_path.is_file():
        raise FileNotFoundError(
            f"Input receptor PDB not found: "
            f"{pdb_path}"
        )

    output_path = (
        work_dir
        / f"{job['job_name']}_reference_ligand.maegz"
    )

    # Reuse an existing extraction unless force=True.
    if (
        output_path.is_file()
        and output_path.stat().st_size > 0
        and not config.get("force", False)
    ):
        logger.info(
            f"Using existing reference ligand: "
            f"{output_path}"
        )

        return output_path

    # Read the receptor/ligand complex.
    entries = list(
        structure.StructureReader(
            str(pdb_path)
        )
    )

    if not entries:
        raise RuntimeError(
            f"No structure found in PDB: "
            f"{pdb_path}"
        )

    st = entries[0]

    # Find atoms belonging to the requested
    # PDB residue name.
    ligand_atoms = [
        atom.index
        for atom in st.atom
        if (atom.pdbres or "").strip().upper()
        == ligand_resname
    ]

    if not ligand_atoms:
        raise RuntimeError(
            f"Reference ligand '{ligand_resname}' "
            f"was not found in {pdb_path.name}."
        )

    n_heavy = heavy_atoms(
        st,
        ligand_atoms,
    )

    if n_heavy == 0:
        raise RuntimeError(
            f"Reference ligand '{ligand_resname}' "
            f"in {pdb_path.name} contains no heavy atoms."
        )

    ligand_st = st.extract(
        ligand_atoms
    )

    ligand_st.title = (
        f"{job['job_name']}_reference_ligand"
    )

    ligand_st.property[
        "s_user_reference_complex"
    ] = job["job_name"]

    ligand_st.property[
        "s_user_reference_resname"
    ] = ligand_resname

    with structure.StructureWriter(
        str(output_path)
    ) as writer:
        writer.append(ligand_st)

    logger.info(
        f"Extracted reference ligand "
        f"'{ligand_resname}' from "
        f"{pdb_path.name} "
        f"({n_heavy} heavy atoms)."
    )

    logger.info(
        f"Reference ligand saved: "
        f"{output_path}"
    )

    return output_path

@register_task(
    "prepare_glideem_receptor",
    category="GlideEM",
    description=(
        "Prepare a receptor using Protein Preparation Wizard."
    ),
)
def prepare_glideem_receptor(
    backend,
    ligand,
    config,
    **kwargs,
):
    job = kwargs["job"]
    work = Path(kwargs["work_dir"])

    source = (
        Path(job["receptor"])
        .expanduser()
        .resolve()
    )

    reference_smiles = (
        job.get("reference_smiles", "")
        or ""
    ).strip()

    output = (
        work /
        f"{job['job_name']}_prepared.maegz"
    )

    if (
        output.is_file()
        and output.stat().st_size > 0
        and not config.get("force", False)
    ):
        logger.info(
            f"Using existing prepared receptor: {output}"
        )
        return output

    prep_input = source

    # --------------------------------------------------
    # SMILES route:
    # remove the reference ligand before PrepWizard.
    #
    # GlideEM will prepare the ligand independently from
    # the authoritative SMILES.
    # --------------------------------------------------

    if reference_smiles:

        ligand_resname = (
            job["reference_ligand"]
            .strip()
            .upper()
        )

        entries = list(
            structure.StructureReader(
                str(source)
            )
        )

        if len(entries) != 1:
            raise RuntimeError(
                f"Expected one structure in {source}, "
                f"found {len(entries)}"
            )

        input_st = entries[0]

        _, ligand_atoms = get_reference_ligand_atoms(
            input_st,
            ligand_resname,
        )

        ligand_atom_set = set(
            ligand_atoms
        )

        receptor_atoms = [
            atom.index
            for atom in input_st.atom
            if atom.index not in ligand_atom_set
        ]

        if not receptor_atoms:
            raise RuntimeError(
                f"No receptor atoms remain after removing "
                f"'{ligand_resname}' from {source}"
            )

        apo_st = input_st.extract(
            receptor_atoms
        )

        prep_input = (
            work /
            f"{job['job_name']}_receptor_without_reference.maegz"
        )

        with structure.StructureWriter(
            str(prep_input)
        ) as writer:
            writer.append(apo_st)

        logger.info(
            f"[{job['job_name']}] Reference SMILES supplied; "
            f"removed '{ligand_resname}' from the receptor "
            f"before PrepWizard"
        )

        logger.info(
            f"[{job['job_name']}] Apo receptor saved: "
            f"{prep_input}"
        )

    # --------------------------------------------------
    # Normal route:
    # retain the reference ligand in the complex.
    # --------------------------------------------------

    else:
        logger.info(
            f"[{job['job_name']}] No reference SMILES supplied; "
            f"preparing the original receptor-ligand complex"
        )

    prep = (
        config
        .get("glideem", {})
        .get("prepare_receptor", {})
    )

    cmd = [
        backend.exe("prepwizard"),
        prep_input,
        output,
        "-WAIT",
        "-HOST",
        backend.host,
        "-fillsidechains",
        "-disulfides",
        "-rehtreat",
        "-propka_pH",
        str(prep.get("pH", 7.4)),
        "-rmsd",
        str(prep.get("rmsd", 0.3)),
    ]

    backend.run(
        cmd,
        work,
        work / f"{job['job_name']}_prepwizard.log",
    )

    if (
        not output.is_file()
        or output.stat().st_size == 0
    ):
        raise RuntimeError(
            f"PrepWizard did not create {output}"
        )

    return output

@register_task(
    "run_glideem",
    category="GlideEM",
    description="Run Schrodinger GlideEM pose placement."
)
def run_glideem(
    backend,
    ligand,
    config,
    **kwargs,
):
    job = kwargs["job"]

    work = Path(
        kwargs["work_dir"]
    )

    receptor = (
        Path(kwargs["prepared_receptor"])
        .expanduser()
        .resolve()
    )

    density_map = (
        Path(job["density_map"])
        .expanduser()
        .resolve()
    )

    if not receptor.is_file():
        raise FileNotFoundError(
            f"Prepared receptor not found: {receptor}"
        )

    if not density_map.is_file():
        raise FileNotFoundError(
            f"Density map not found: {density_map}"
        )

    # The manifest contains the PDB residue name,
    # for example LIG.
    ligand_resname = (
        job["reference_ligand"]
        .strip()
        .upper()
    )

    if not ligand_resname:
        raise ValueError(
            f"Empty reference_ligand value for "
            f"{job['job_name']}"
        )

    # GlideEM help uses this exact ASL style:
    #
    #     res.ptype 8PR
    #
    # so use the manifest residue name in the same way.
    reference_smiles = (
        job.get("reference_smiles", "")
        or ""
    ).strip()

    ligand_asl = None
    docking_site_xyz = None

    if reference_smiles:

        source_structure = (
            Path(job["receptor"])
            .expanduser()
            .resolve()
        )

        docking_site_xyz = get_reference_ligand_center(
            source_structure,
            ligand_resname,
        )

        logger.info(
            f"[{job['job_name']}] Using "
            f"reference SMILES for GlideEM ligand preparation"
        )

        logger.debug(
            f"[{job['job_name']}] Experimental docking center: "
            f"{docking_site_xyz[0]:.3f}, "
            f"{docking_site_xyz[1]:.3f}, "
            f"{docking_site_xyz[2]:.3f}"
        )


    else:

        ligand_asl = (
            f"res.ptype {ligand_resname}"
        )

        logger.info(
            f"[{job['job_name']}] Using prepared complex "
            f"ligand ASL: '{ligand_asl}'"
        )

    glideem_cfg = config.get(
        "glideem",
        {}
    )

    mode = glideem_cfg.get(
        "mode",
        "exhaustive",
    )

    nposes = int(
        glideem_cfg.get(
            "nposes",
            5,
        )
    )

    job_name = (
        f"{job['job_name']}_glideem"
    )

    ligand_source = (
        "SMILES"
        if reference_smiles
        else f"ASL '{ligand_asl}'"
    )

    logger.info(
        f"Starting GlideEM for {job['job_name']}: "
        f"ligand_source={ligand_source}, "
        f"mode={mode}, nposes={nposes}"
    )

    backend.run_glideem(
        receptor=receptor,
        density_map=density_map,
        job_name=job_name,
        cwd=work,
        mode=mode,
        nposes=nposes,
        ligand_asl=ligand_asl,
        ligand_smiles=(
            reference_smiles or None
        ),
        docking_site_xyz=docking_site_xyz,
        host=glideem_cfg.get(
            "host",
            backend.host,
        ),
    )

    job_log = (
        work /
        f"{job_name}.log"
    )

    if not job_log.is_file():
        raise RuntimeError(
            f"GlideEM Job Control log was not created: "
            f"{job_log}"
        )

    job_log_text = job_log.read_text(
        errors="replace"
    )

    failure_markers = [
        "Traceback (most recent call last)",
        "RuntimeError:",
        "ERROR:",
    ]

    detected_failures = [
        marker
        for marker in failure_markers
        if marker in job_log_text
    ]

    if detected_failures:

        tail = "\n".join(
            job_log_text.splitlines()[-80:]
        )

        raise RuntimeError(
            f"GlideEM job '{job_name}' failed. "
            f"Detected markers: "
            f"{', '.join(detected_failures)}\n"
            f"Job log: {job_log}\n"
            f"Last 80 lines:\n{tail}"
        )
    output_file = (
        work /
        f"{job_name}-out.maegz"
    )

    if not output_file.is_file():
        raise RuntimeError(
            f"GlideEM completed without producing the expected "
            f"pose file: {output_file}"
        )

    if output_file.stat().st_size == 0:
        raise RuntimeError(
            f"GlideEM pose file is empty: {output_file}"
        )

    logger.info(
        f"GlideEM completed successfully for "
        f"{job['job_name']}: {output_file}"
    )

    return output_file


@register_task(
    "analyse_glideem_poses",
    category="GlideEM",
    description=(
        "Analyse GlideEM pose scores, RMSD distribution, "
        "pose clustering and PyMOL visualization."
    ),
)
def analyse_glideem_poses(
    backend,
    ligand,
    config,
    **kwargs,
):
    archive = (
        Path(kwargs["archive"])
        .expanduser()
        .resolve()
    )

    glideem_csv = (
        Path(kwargs["glideem_csv"])
        .expanduser()
        .resolve()
    )

    job = kwargs["job"]

    output_dir = (
        Path(kwargs["output_dir"])
        .expanduser()
        .resolve()
    )

    cryoem_reference = (
        Path(kwargs["cryoem_reference"])
        .expanduser()
        .resolve()
    )

    ligand_resname = (
        job["reference_ligand"]
        .strip()
        .upper()
    )

    prepared_receptor = (
        Path(kwargs["prepared_receptor"])
        .expanduser()
        .resolve()
    )

    density_map = (
        Path(job["density_map"])
        .expanduser()
        .resolve()
    )


    reference_smiles = (
        job.get("reference_smiles", "")
        or ""
    ).strip()

    ligand_from_smiles = bool(
        reference_smiles
    )

    analysis_cfg = config.get(
        "glideem_analysis",
        {},
    )

    cluster_cutoff = float(
        analysis_cfg.get(
            "cluster_rmsd_cutoff",
            2.0,
        )
    )

    export_pymol = bool(
        analysis_cfg.get(
            "export_pymol",
            True,
        )
    )

    if not archive.is_file():
        raise FileNotFoundError(
            f"GlideEM archive not found: {archive}"
        )

    if not glideem_csv.is_file():
        raise FileNotFoundError(
            f"GlideEM metrics CSV not found: {glideem_csv}"
        )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    # --------------------------------------------------
    # Read the GlideEM score table
    # --------------------------------------------------

    metrics = pd.read_csv(
        glideem_csv
    )

    required_columns = {
        "Pose",
        "Glide Density Score",
        "Glide Docking Score",
        "Ligand RMSD",
    }

    missing_columns = (
        required_columns
        - set(metrics.columns)
    )

    if missing_columns:
        raise RuntimeError(
            f"GlideEM CSV is missing columns: "
            f"{sorted(missing_columns)}"
        )

    metrics = metrics.copy()

    # CSV poses are zero-based; pipeline/PyMOL names are one-based.
    metrics["pose_number"] = (
        metrics["Pose"].astype(int)
        + 1
    )

    # --------------------------------------------------
    # Read ligand coordinates from each GlideEM CT
    # --------------------------------------------------

    pose_structures = []
    pose_coordinates = []
    expected_elements = None

    for ct_number, st in enumerate(
        structure.StructureReader(
            str(archive)
        ),
        start=1,
    ):
        ligand_indices = find_glideem_ligand_indices(
            st,
            ligand_resname=ligand_resname,
            ligand_from_smiles=ligand_from_smiles,
        )

        if not ligand_indices:
            logger.warning(
                "[%s] GlideEM CT %d contains no usable ligand; "
                "excluded from ensemble analysis",
                job["job_name"],
                ct_number,
            )
            continue

        coordinates, elements = (
            ligand_heavy_atom_coordinates(
                st,
                ligand_indices,
            )
        )

        if expected_elements is None:
            expected_elements = elements

        elif elements != expected_elements:
            raise RuntimeError(
                f"GlideEM pose {ct_number} has a different "
                "heavy-atom element sequence. Symmetry-aware "
                "atom mapping is required before RMSD comparison."
            )

        ligand_st = st.extract(
            ligand_indices
        )

        pose_structures.append(
            {
                "pose_number": ct_number,
                "structure": ligand_st,
            }
        )

        pose_coordinates.append(
            coordinates
        )

    nposes = len(
        pose_coordinates
    )

    if nposes == 0:
        raise RuntimeError(
            f"No usable GlideEM poses found in {archive}"
        )

    # Keep only score rows corresponding to extracted poses.
    extracted_pose_numbers = {
        record["pose_number"]
        for record in pose_structures
    }

    metrics = metrics[
        metrics["pose_number"].isin(
            extracted_pose_numbers
        )
    ].copy()

    # --------------------------------------------------
    # All-vs-all RMSD
    # --------------------------------------------------

    rmsd_matrix = np.zeros(
        (nposes, nposes),
        dtype=float,
    )

    for i in range(nposes):

        for j in range(i + 1, nposes):

            value = coordinate_frame_rmsd(
                pose_coordinates[i],
                pose_coordinates[j],
            )

            rmsd_matrix[i, j] = value
            rmsd_matrix[j, i] = value

    pose_labels = [
        f"pose{record['pose_number']:02d}"
        for record in pose_structures
    ]

    rmsd_df = pd.DataFrame(
        rmsd_matrix,
        index=pose_labels,
        columns=pose_labels,
    )

    rmsd_df.to_csv(
        output_dir /
        "pose_pairwise_rmsd.csv"
    )

    # --------------------------------------------------
    # Hierarchical clustering
    # --------------------------------------------------

    if nposes > 1:

        condensed_distances = squareform(
            rmsd_matrix,
            checks=False,
        )

        linkage_matrix = linkage(
            condensed_distances,
            method="average",
        )

        cluster_ids = fcluster(
            linkage_matrix,
            t=cluster_cutoff,
            criterion="distance",
        )

    else:
        linkage_matrix = None
        cluster_ids = np.ones(
            1,
            dtype=int,
        )

    cluster_by_pose = {
        record["pose_number"]: int(cluster_id)
        for record, cluster_id in zip(
            pose_structures,
            cluster_ids,
        )
    }

    metrics["cluster_id"] = (
        metrics["pose_number"]
        .map(cluster_by_pose)
    )

    metrics.to_csv(
        output_dir /
        "glideem_pose_metrics.csv",
        index=False,
    )

    cluster_summary = (
        metrics
        .groupby("cluster_id", dropna=False)
        .agg(
            n_poses=("pose_number", "size"),
            mean_density_score=(
                "Glide Density Score",
                "mean",
            ),
            best_density_score=(
                "Glide Density Score",
                "min",
            ),
            mean_docking_score=(
                "Glide Docking Score",
                "mean",
            ),
            best_docking_score=(
                "Glide Docking Score",
                "min",
            ),
            mean_reference_rmsd=(
                "Ligand RMSD",
                "mean",
            ),
            min_reference_rmsd=(
                "Ligand RMSD",
                "min",
            ),
            max_reference_rmsd=(
                "Ligand RMSD",
                "max",
            ),
        )
        .reset_index()
    )

    cluster_summary.to_csv(
        output_dir /
        "pose_clusters.csv",
        index=False,
    )

    # --------------------------------------------------
    # Plots
    # --------------------------------------------------

    save_clustered_score_vs_rmsd(
        metrics=metrics,
        score_column="Glide Density Score",
        output_file=(
            output_dir /
            "density_score_vs_rmsd.png"
        ),
        title=(
            f"{job['job_name']}: "
            "density score vs RMSD"
        ),
        y_label="Glide Density Score",
    )

    save_clustered_score_vs_rmsd(
        metrics=metrics,
        score_column="Glide Docking Score",
        output_file=(
            output_dir /
            "docking_score_vs_rmsd.png"
        ),
        title=(
            f"{job['job_name']}: "
            "docking score vs RMSD"
        ),
        y_label="Glide Docking Score",
    )

    save_rmsd_heatmap(
        rmsd_matrix=rmsd_matrix,
        pose_labels=pose_labels,
        output_file=(
            output_dir /
            "pose_rmsd_heatmap.png"
        ),
        title=(
            f"{job['job_name']}: "
            "pose-to-pose RMSD"
        ),
    )

    if linkage_matrix is not None:
        save_cluster_dendrogram(
            linkage_matrix=linkage_matrix,
            pose_labels=pose_labels,
            cluster_cutoff=cluster_cutoff,
            output_file=(
                output_dir /
                "pose_cluster_dendrogram.png"
            ),
            title=(
                f"{job['job_name']}: "
                "GlideEM pose clustering"
            ),
        )

    # --------------------------------------------------
    # PyMOL export
    # --------------------------------------------------

    pymol_script = None

    if export_pymol:

        pymol_script = export_glideem_pymol_ensemble(
            pose_structures=pose_structures,
            cluster_by_pose=cluster_by_pose,
            cryoem_reference=cryoem_reference,
            receptor=prepared_receptor,
            density_map=density_map,
            output_dir=output_dir / "pymol",
            job_name=job["job_name"],
        )

    logger.info(
        "[%s] Analysed %d GlideEM poses in %d clusters",
        job["job_name"],
        nposes,
        len(set(cluster_ids)),
    )

    return {
        "job_name": job["job_name"],
        "requested_nposes": int(
            config.get(
                "glideem",
                {},
            ).get(
                "nposes",
                5,
            )
        ),
        "cluster_rmsd_cutoff": cluster_cutoff,
        "metrics": (
            output_dir /
            "glideem_pose_metrics.csv"
        ),
        "clusters": (
            output_dir /
            "pose_clusters.csv"
        ),
        "pairwise_rmsd": (
            output_dir /
            "pose_pairwise_rmsd.csv"
        ),
        "pymol_script": pymol_script,
    }

@register_task(
    "summarise_glideem_ensembles",
    category="GlideEM",
    description=(
        "Summarise GlideEM pose-space ambiguity "
        "across multiple structures and density maps."
    ),
)
def summarise_glideem_ensembles(
    backend,
    ligand,
    config,
    **kwargs,
):
    analysis_results = list(
        kwargs["analysis_results"]
    )

    output_dir = (
        Path(kwargs["output_dir"])
        .expanduser()
        .resolve()
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not analysis_results:
        raise RuntimeError(
            "No successful GlideEM ensemble analyses "
            "were supplied for campaign summarisation."
        )

    summary_rows = []
    all_pose_metrics = []
    cluster_composition_rows = []

    # --------------------------------------------------
    # Read per-structure outputs
    # --------------------------------------------------

    for result in analysis_results:

        job_name = result["job_name"]

        metrics_file = Path(
            result["metrics"]
        )

        clusters_file = Path(
            result["clusters"]
        )

        pairwise_file = Path(
            result["pairwise_rmsd"]
        )

        if not metrics_file.is_file():
            logger.warning(
                "[%s] GlideEM metrics file missing; "
                "excluded from campaign summary: %s",
                job_name,
                metrics_file,
            )
            continue

        if not clusters_file.is_file():
            logger.warning(
                "[%s] GlideEM cluster file missing; "
                "excluded from campaign summary: %s",
                job_name,
                clusters_file,
            )
            continue

        if not pairwise_file.is_file():
            logger.warning(
                "[%s] Pairwise RMSD file missing; "
                "excluded from campaign summary: %s",
                job_name,
                pairwise_file,
            )
            continue

        metrics = pd.read_csv(
            metrics_file
        )

        clusters = pd.read_csv(
            clusters_file
        )

        pairwise = pd.read_csv(
            pairwise_file,
            index_col=0,
        )

        if metrics.empty:
            logger.warning(
                "[%s] GlideEM metrics table is empty; "
                "excluded from campaign summary",
                job_name,
            )
            continue

        requested_nposes = int(
            result["requested_nposes"]
        )

        returned_nposes = len(
            metrics
        )

        cluster_cutoff = float(
            result["cluster_rmsd_cutoff"]
        )

        cluster_counts = (
            metrics["cluster_id"]
            .dropna()
            .astype(int)
            .value_counts()
            .sort_index()
        )

        n_clusters = len(
            cluster_counts
        )

        if returned_nposes:
            cluster_fractions = (
                cluster_counts
                / returned_nposes
            )
        else:
            cluster_fractions = pd.Series(
                dtype=float
            )

        dominant_cluster_fraction = (
            float(cluster_fractions.max())
            if not cluster_fractions.empty
            else np.nan
        )

        # Normalised Shannon entropy:
        #
        # 0 = one occupied family
        # 1 = occupancy evenly distributed among families
        if n_clusters <= 1:
            cluster_entropy = 0.0

        else:
            fractions = (
                cluster_fractions
                .to_numpy(float)
            )

            raw_entropy = float(
                -np.sum(
                    fractions
                    * np.log(fractions)
                )
            )

            cluster_entropy = (
                raw_entropy
                / math.log(n_clusters)
            )

        # Pairwise upper triangle, excluding diagonal.
        pairwise_values = (
            pairwise
            .to_numpy(dtype=float)
        )

        if pairwise_values.shape[0] > 1:
            upper_indices = np.triu_indices(
                pairwise_values.shape[0],
                k=1,
            )

            pairwise_distances = (
                pairwise_values[
                    upper_indices
                ]
            )

            pairwise_distances = (
                pairwise_distances[
                    np.isfinite(
                        pairwise_distances
                    )
                ]
            )

        else:
            pairwise_distances = np.array(
                [],
                dtype=float,
            )

        median_pairwise_rmsd = (
            float(
                np.median(
                    pairwise_distances
                )
            )
            if pairwise_distances.size
            else 0.0
        )

        p90_pairwise_rmsd = (
            float(
                np.percentile(
                    pairwise_distances,
                    90,
                )
            )
            if pairwise_distances.size
            else 0.0
        )

        reference_rmsd = (
            metrics["Ligand RMSD"]
            .to_numpy(float)
        )

        density_scores = (
            metrics["Glide Density Score"]
            .to_numpy(float)
        )

        docking_scores = (
            metrics["Glide Docking Score"]
            .to_numpy(float)
        )

        if (
            len(metrics) >= 2
            and np.ptp(reference_rmsd) > 0
            and np.ptp(density_scores) > 0
        ):
            density_rmsd_result = spearmanr(
                reference_rmsd,
                density_scores,
            )

            density_rmsd_spearman = float(
                density_rmsd_result.statistic
            )

            density_rmsd_p = float(
                density_rmsd_result.pvalue
            )

        else:
            density_rmsd_spearman = np.nan
            density_rmsd_p = np.nan

        summary_rows.append(
            {
                "source_name": job_name,
                "requested_nposes": requested_nposes,
                "returned_nposes": returned_nposes,
                "return_fraction": (
                    returned_nposes
                    / requested_nposes
                    if requested_nposes
                    else np.nan
                ),
                "cluster_rmsd_cutoff": cluster_cutoff,
                "n_clusters": n_clusters,
                "dominant_cluster_fraction": (
                    dominant_cluster_fraction
                ),
                "cluster_entropy": cluster_entropy,
                "median_reference_rmsd": float(
                    np.median(
                        reference_rmsd
                    )
                ),
                "reference_rmsd_iqr": float(
                    np.percentile(
                        reference_rmsd,
                        75,
                    )
                    - np.percentile(
                        reference_rmsd,
                        25,
                    )
                ),
                "max_reference_rmsd": float(
                    np.max(
                        reference_rmsd
                    )
                ),
                "median_pairwise_rmsd": (
                    median_pairwise_rmsd
                ),
                "p90_pairwise_rmsd": (
                    p90_pairwise_rmsd
                ),
                "density_score_range": float(
                    np.max(
                        density_scores
                    )
                    - np.min(
                        density_scores
                    )
                ),
                "docking_score_range": float(
                    np.max(
                        docking_scores
                    )
                    - np.min(
                        docking_scores
                    )
                ),
                "density_rmsd_spearman": (
                    density_rmsd_spearman
                ),
                "density_rmsd_p": (
                    density_rmsd_p
                ),
            }
        )

        metrics_for_campaign = (
            metrics.copy()
        )

        metrics_for_campaign[
            "source_name"
        ] = job_name

        all_pose_metrics.append(
            metrics_for_campaign
        )

        for cluster_id, count in (
            cluster_counts.items()
        ):
            cluster_composition_rows.append(
                {
                    "source_name": job_name,
                    "cluster_id": int(
                        cluster_id
                    ),
                    "n_poses": int(
                        count
                    ),
                    "cluster_fraction": float(
                        count
                        / returned_nposes
                    ),
                }
            )

    if not summary_rows:
        raise RuntimeError(
            "No valid GlideEM ensemble analyses "
            "were available for campaign summarisation."
        )

    summary_df = pd.DataFrame(
        summary_rows
    )

    summary_df = summary_df.sort_values(
        [
            "cluster_entropy",
            "median_pairwise_rmsd",
        ],
        ascending=[
            False,
            False,
        ],
    ).reset_index(
        drop=True
    )

    summary_df.to_csv(
        output_dir /
        "glideem_ensemble_summary.csv",
        index=False,
    )

    all_metrics_df = pd.concat(
        all_pose_metrics,
        ignore_index=True,
    )

    all_metrics_df.to_csv(
        output_dir /
        "all_glideem_pose_metrics.csv",
        index=False,
    )

    composition_df = pd.DataFrame(
        cluster_composition_rows
    )

    composition_df.to_csv(
        output_dir /
        "cluster_composition.csv",
        index=False,
    )

    # --------------------------------------------------
    # Campaign plots
    # --------------------------------------------------

    job_order = (
        summary_df[
            "source_name"
        ]
        .tolist()
    )

    save_pose_ambiguity_map(
        summary_df=summary_df,
        output_file=(
            output_dir /
            "pose_ambiguity_map.png"
        ),
    )

    save_reference_rmsd_distributions(
        all_metrics_df=all_metrics_df,
        job_order=job_order,
        output_file=(
            output_dir /
            "rmsd_distribution_by_structure.png"
        ),
    )

    save_cluster_composition(
        composition_df=composition_df,
        summary_df=summary_df,
        job_order=job_order,
        output_file=(
            output_dir /
            "cluster_composition_by_structure.png"
        ),
    )

    save_density_rmsd_small_multiples(
        all_metrics_df=all_metrics_df,
        job_order=job_order,
        output_file=(
            output_dir /
            "density_rmsd_small_multiples.png"
        ),
        top_n_clusters=10,
    )

    logger.info(
        "Generated GlideEM ensemble summary for "
        "%d structures",
        len(summary_df),
    )

    return summary_df

@register_task("build_glide_grid", category="Glide", description="Build a Glide grid around a GlideEM reference pose")
def build_glide_grid(backend, ligand, config, **kwargs):
    job = kwargs["complex_job"]
    grid_cfg = config["grid"]
    if job["grid"].exists() and not config.get("force", False):
        return job["grid"]
    x, y, z = job["center"]
    inp = job["dir"] / f"{job['name']}_grid.in"
    inp.write_text(
        f"GRIDFILE {job['grid']}\nRECEP_FILE {job['receptor']}\n"
        f"GRID_CENTER {x:.4f},{y:.4f},{z:.4f}\n"
        f"INNERBOX {grid_cfg.get('innerbox', 10)},{grid_cfg.get('innerbox', 10)},{grid_cfg.get('innerbox', 10)}\n"
        f"OUTERBOX {grid_cfg.get('outerbox', 30)},{grid_cfg.get('outerbox', 30)},{grid_cfg.get('outerbox', 30)}\n"
    )
    backend.run_glide(
        inp,
        job["dir"],
        job["dir"] / f"{job['name']}_grid.log",
        host=grid_cfg.get(
            "host",
            backend.host,
        ),
    )
    if not job["grid"].exists():
        raise RuntimeError(f"Grid missing: {job['grid']}")
    return job["grid"]


@register_task(
    "dock_mcss",
    category="Glide",
    description="Core-constrained MCSS docking using GlideEM pose",
)
def dock_mcss(
    backend,
    ligand,
    config,
    **kwargs,
):
    job = kwargs["complex_job"]

    ligand_file = (
        Path(kwargs["ligand_file"])
        .expanduser()
        .resolve()
    )

    out_dir = (
        Path(kwargs["mcss_root"])
        / job["name"]
    )

    out_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    inp = (
        out_dir /
        f"{job['name']}_mcss.in"
    )

    pv = (
        out_dir /
        f"{job['name']}_mcss_pv.maegz"
    )

    launcher_log = (
        out_dir /
        f"{job['name']}_mcss_launcher.log"
    )

    mcss = config.get(
        "mcss",
        {}
    )

    # --------------------------------------------------
    # Validate inputs
    # --------------------------------------------------

    grid_file = (
        Path(job["grid"])
        .expanduser()
        .resolve()
    )

    if not grid_file.is_file():
        raise RuntimeError(
            f"Cannot start MCSS docking for "
            f"{job['name']}: Glide grid does not exist: "
            f"{grid_file}"
        )

    if grid_file.stat().st_size == 0:
        raise RuntimeError(
            f"Cannot start MCSS docking for "
            f"{job['name']}: Glide grid is empty: "
            f"{grid_file}"
        )

    if not ligand_file.is_file():
        raise FileNotFoundError(
            f"Prepared ligand library not found: "
            f"{ligand_file}"
        )

    reference = (
        Path(job["reference"])
        .expanduser()
        .resolve()
    )

    if not reference.is_file():
        raise FileNotFoundError(
            f"Reference ligand not found: "
            f"{reference}"
        )

    # --------------------------------------------------
    # Write Glide docking input
    # --------------------------------------------------

    inp.write_text(
        "CORECONS_FALLBACK True\n"
        "CORE_DEFINITION mcssmarts\n"
        "CORE_FILTER True\n"
        "CORE_RESTRAIN True\n"
        f"GRIDFILE {grid_file}\n"
        f"LIGANDFILE {ligand_file}\n"
        f"NREPORT {mcss.get('nreport', 10000)}\n"
        f"POSES_PER_LIG "
        f"{mcss.get('poses_per_ligand', 1)}\n"
        f"POSTDOCK_NPOSE "
        f"{mcss.get('postdock_npose', 1)}\n"
        f"PRECISION {mcss.get('precision', 'SP')}\n"
        f"REF_LIGAND_FILE {reference}\n"
        "USE_REF_LIGAND True\n"
    )

    # --------------------------------------------------
    # Run Glide
    # --------------------------------------------------

    if (
        not pv.exists()
        or config.get("force", False)
    ):
        logger.info(
            f"[{job['name']}] Launching "
            f"core-constrained Glide docking"
        )

        backend.run_glide(
            inp,
            out_dir,
            launcher_log,
            host=mcss.get(
                "host",
                backend.host,
            ),
        )

    if not pv.exists():
        raise RuntimeError(
            f"MCSS PoseViewer missing after Glide "
            f"docking: {pv}"
        )

    return pv

@register_task(
    "export_mcss_pymol",
    category="GlideEM",
    description=(
        "Export top-scoring MCSS poses "
        "for PyMOL visualization."
    ),
)
def export_mcss_pymol(
    backend,
    ligand,
    config,
    **kwargs,
):
    complex_job = kwargs["complex_job"]

    poseviewer = (
        Path(kwargs["poseviewer"])
        .expanduser()
        .resolve()
    )

    top_n = int(
        kwargs.get(
            "top_n",
            5,
        )
    )

    cryoem_reference = (
        Path(kwargs["cryoem_reference"])
        .expanduser()
        .resolve()
    )

    glideem_reference = (
        Path(complex_job["reference"])
        .expanduser()
        .resolve()
    )

    pymol_cfg = config.get(
        "pymol",
        {},
    )

    binding_site_cutoff = float(
        pymol_cfg.get(
            "binding_site_cutoff",
            5.0,
        )
    )

    if not poseviewer.is_file():
        raise FileNotFoundError(
            f"MCSS PoseViewer not found: {poseviewer}"
        )

    if not cryoem_reference.is_file():
        raise FileNotFoundError(
            "Original cryo-EM reference ligand "
            f"not found: {cryoem_reference}"
        )

    if not glideem_reference.is_file():
        raise FileNotFoundError(
            "GlideEM reference pose not found: "
            f"{glideem_reference}"
        )

    if top_n < 1:
        raise ValueError(
            "PyMOL top_n must be at least 1."
        )

    name = complex_job["name"]

    out_dir = (
        poseviewer.parent
        / "pymol"
    )

    receptor = None
    poses = []

    for st in structure.StructureReader(
        str(poseviewer)
    ):
        if st.property.get(
            "b_glide_receptor",
            False,
        ):
            receptor = st
            continue

        score = glide_score(
            st
        )

        if score is None:
            continue

        cid = compound_id(
            st
        )

        poses.append(
            (
                float(score),
                cid,
                st,
            )
        )

    if receptor is None:
        raise RuntimeError(
            f"No receptor CT found in {poseviewer}"
        )

    if not poses:
        raise RuntimeError(
            "No scored docked poses found in "
            f"{poseviewer}"
        )

    # More-negative GlideScore is better.
    poses.sort(
        key=lambda item: item[0]
    )

    top_poses = poses[
        :top_n
    ]

    ranked_poses = [
        {
            "rank": rank,
            "score": score,
            "compound_id": cid,
            "structure": st,
        }
        for rank, (
            score,
            cid,
            st,
        ) in enumerate(
            top_poses,
            start=1,
        )
    ]

    script_file = export_mcss_pose_bundle(
        receptor=receptor,
        ranked_poses=ranked_poses,
        cryoem_reference=cryoem_reference,
        glideem_reference=glideem_reference,
        output_dir=out_dir,
        binding_site_cutoff=(
            binding_site_cutoff
        ),
    )

    logger.debug(
        "[%s] Exported top %d MCSS poses "
        "for PyMOL to %s",
        name,
        len(ranked_poses),
        out_dir,
    )

    return script_file