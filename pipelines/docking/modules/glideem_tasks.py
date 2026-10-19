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

    output_root = Path(config["output_dir"]).expanduser().resolve()

    output_dir = output_root / "ligprep"
    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    output_file = (
        output_dir /
        "activity_ligands_ligprep.maegz"
    )

    # Reuse existing output unless explicitly forcing rerun
    if (
        output_file.is_file()
        and output_file.stat().st_size > 0
        and not config.get("force", False)
    ):
        logger.info(
            f"Using existing LigPrep output: {output_file}"
        )

        return output_file

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

    cmd = [
        backend.exe("ligprep"),
        "-isd",
        str(input_path),
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
    # Density score versus reference RMSD
    # --------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(7, 6)
    )

    scatter = ax.scatter(
        metrics["Ligand RMSD"],
        metrics["Glide Density Score"],
        c=metrics["cluster_id"],
        alpha=0.85,
    )

    ax.set_xlabel(
        "RMSD from input cryo-EM ligand pose (Å)"
    )

    ax.set_ylabel(
        "Glide Density Score"
    )

    ax.set_title(
        f"{job['job_name']}: density score vs RMSD"
    )

    fig.colorbar(
        scatter,
        ax=ax,
        label="Pose cluster",
    )

    fig.tight_layout()

    fig.savefig(
        output_dir /
        "density_score_vs_rmsd.png",
        dpi=200,
    )

    plt.close(fig)

    # --------------------------------------------------
    # Docking score versus reference RMSD
    # --------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(7, 6)
    )

    scatter = ax.scatter(
        metrics["Ligand RMSD"],
        metrics["Glide Docking Score"],
        c=metrics["cluster_id"],
        alpha=0.85,
    )

    ax.set_xlabel(
        "RMSD from input cryo-EM ligand pose (Å)"
    )

    ax.set_ylabel(
        "Glide Docking Score"
    )

    ax.set_title(
        f"{job['job_name']}: docking score vs RMSD"
    )

    fig.colorbar(
        scatter,
        ax=ax,
        label="Pose cluster",
    )

    fig.tight_layout()

    fig.savefig(
        output_dir /
        "docking_score_vs_rmsd.png",
        dpi=200,
    )

    plt.close(fig)

    # --------------------------------------------------
    # RMSD heatmap
    # --------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(8, 7)
    )

    image = ax.imshow(
        rmsd_matrix,
        interpolation="nearest",
        aspect="auto",
    )

    ax.set_xticks(
        np.arange(nposes)
    )

    ax.set_yticks(
        np.arange(nposes)
    )

    ax.set_xticklabels(
        pose_labels,
        rotation=90,
    )

    ax.set_yticklabels(
        pose_labels,
    )

    ax.set_title(
        f"{job['job_name']}: pose-to-pose RMSD"
    )

    fig.colorbar(
        image,
        ax=ax,
        label="Heavy-atom RMSD (Å)",
    )

    fig.tight_layout()

    fig.savefig(
        output_dir /
        "pose_rmsd_heatmap.png",
        dpi=200,
    )

    plt.close(fig)

    # --------------------------------------------------
    # Dendrogram
    # --------------------------------------------------

    if linkage_matrix is not None:

        fig, ax = plt.subplots(
            figsize=(9, 5)
        )

        dendrogram(
            linkage_matrix,
            labels=pose_labels,
            ax=ax,
        )

        ax.axhline(
            cluster_cutoff,
            linestyle="--",
            label=(
                f"Cluster cutoff "
                f"{cluster_cutoff:.2f} Å"
            ),
        )

        ax.set_ylabel(
            "Pose RMSD (Å)"
        )

        ax.set_title(
            f"{job['job_name']}: GlideEM pose clustering"
        )

        ax.legend()

        fig.tight_layout()

        fig.savefig(
            output_dir /
            "pose_cluster_dendrogram.png",
            dpi=200,
        )

        plt.close(fig)

    # --------------------------------------------------
    # PyMOL export
    # --------------------------------------------------

    pymol_script = None

    if export_pymol:

        pymol_script = export_glideem_pymol_ensemble(
            pose_structures=pose_structures,
            cluster_by_pose=cluster_by_pose,
            cryoem_reference=cryoem_reference,
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
    # Plot 1: pose ambiguity map
    # --------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(9, 7)
    )

    point_sizes = (
        60
        + 18
        * summary_df[
            "returned_nposes"
        ].to_numpy(float)
    )

    scatter = ax.scatter(
        summary_df[
            "median_pairwise_rmsd"
        ],
        summary_df[
            "cluster_entropy"
        ],
        s=point_sizes,
        c=summary_df[
            "dominant_cluster_fraction"
        ],
        alpha=0.85,
        vmin=0.0,
        vmax=1.0,
    )

    for _, row in summary_df.iterrows():

        ax.annotate(
            row["source_name"],
            (
                row[
                    "median_pairwise_rmsd"
                ],
                row[
                    "cluster_entropy"
                ],
            ),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
        )

    ax.set_xlabel(
        "Median pairwise pose RMSD (Å)"
    )

    ax.set_ylabel(
        "Normalised cluster entropy"
    )

    ax.set_title(
        "GlideEM pose-space ambiguity"
    )

    colour_bar = fig.colorbar(
        scatter,
        ax=ax,
    )

    colour_bar.set_label(
        "Dominant-cluster fraction"
    )

    ax.text(
        0.02,
        0.02,
        "Point size = returned GlideEM poses",
        transform=ax.transAxes,
        fontsize=9,
        verticalalignment="bottom",
    )

    fig.tight_layout()

    fig.savefig(
        output_dir /
        "pose_ambiguity_map.png",
        dpi=200,
    )

    plt.close(fig)

    # --------------------------------------------------
    # Plot 2: reference RMSD distributions
    # --------------------------------------------------

    job_order = (
        summary_df["source_name"]
        .tolist()
    )

    rmsd_distributions = [
        all_metrics_df.loc[
            all_metrics_df[
                "source_name"
            ] == job_name,
            "Ligand RMSD",
        ].to_numpy(float)
        for job_name in job_order
    ]

    fig_height = max(
        5,
        0.6 * len(job_order) + 2,
    )

    fig, ax = plt.subplots(
        figsize=(9, fig_height)
    )

    positions = np.arange(
        1,
        len(job_order) + 1,
    )

    violin_parts = ax.violinplot(
        rmsd_distributions,
        positions=positions,
        vert=False,
        showmeans=False,
        showmedians=True,
        showextrema=True,
    )

    for body in violin_parts[
        "bodies"
    ]:
        body.set_alpha(
            0.65
        )

    ax.set_yticks(
        positions
    )

    ax.set_yticklabels(
        job_order
    )

    ax.set_xlabel(
        "RMSD from input cryo-EM ligand pose (Å)"
    )

    ax.set_ylabel(
        "GlideEM job"
    )

    ax.set_title(
        "GlideEM reference-pose RMSD distributions"
    )

    fig.tight_layout()

    fig.savefig(
        output_dir /
        "rmsd_distribution_by_structure.png",
        dpi=200,
    )

    plt.close(fig)

    # --------------------------------------------------
    # Plot 3: cluster composition
    # --------------------------------------------------

    cluster_ids = sorted(
        composition_df[
            "cluster_id"
        ].unique()
    )

    composition_pivot = (
        composition_df
        .pivot(
            index="source_name",
            columns="cluster_id",
            values="cluster_fraction",
        )
        .fillna(0.0)
        .reindex(job_order)
    )

    clusters_per_job = {
        row["source_name"]: int(row["n_clusters"])
        for _, row in summary_df.iterrows()
    }

    ylabels = [
        f"{job_name}  (k={clusters_per_job[job_name]})"
        for job_name in composition_pivot.index
    ]

    fig_height = max(
        5,
        0.55 * len(job_order) + 2,
    )

    fig, ax = plt.subplots(
        figsize=(10, fig_height)
    )

    left = np.zeros(
        len(composition_pivot),
        dtype=float,
    )

    for cluster_id in cluster_ids:

        if cluster_id in composition_pivot.columns:
            values = (
                composition_pivot[
                    cluster_id
                ].to_numpy(float)
            )
        else:
            values = np.zeros(
                len(composition_pivot),
                dtype=float,
            )

        bars = ax.barh(
            composition_pivot.index,
            values,
            left=left,
        )

        # Label sufficiently large segments directly.
        for bar, value, start in zip(
            bars,
            values,
            left,
        ):
            if value >= 0.07:
                ax.text(
                    start + value / 2,
                    bar.get_y() + bar.get_height() / 2,
                    f"C{cluster_id}",
                    ha="center",
                    va="center",
                    fontsize=8,
                )

        left += values

    ax.set_yticks(
        np.arange(
            len(composition_pivot)
        )
    )

    ax.set_yticklabels(
        ylabels
    )

    ax.set_xlim(
        0.0,
        1.0,
    )

    ax.set_xlabel(
        "Fraction of returned GlideEM poses"
    )

    ax.set_ylabel(
        "GlideEM job"
    )

    ax.set_title(
        "GlideEM pose-cluster composition"
    )

    fig.tight_layout()

    fig.savefig(
        output_dir /
        "cluster_composition_by_structure.png",
        dpi=200,
    )

    plt.close(fig)

    # --------------------------------------------------
    # Plot 4: density score vs RMSD small multiples
    # --------------------------------------------------

    n_jobs = len(
        job_order
    )

    n_columns = min(
        2,
        n_jobs,
    )

    n_rows = int(
        math.ceil(
            n_jobs / n_columns
        )
    )

    fig, axes = plt.subplots(
        n_rows,
        n_columns,
        figsize=(
            7 * n_columns,
            5.5 * n_rows,
        ),
        squeeze=False,
    )

    for axis, job_name in zip(
        axes.flat,
        job_order,
    ):

        job_metrics = (
            all_metrics_df[
                all_metrics_df[
                    "source_name"
                ] == job_name
            ]
            .copy()
        )

        # ---------------------------------------------
        # Keep the 10 most populous clusters.
        # Everything else becomes "Other".
        # ---------------------------------------------

        cluster_counts = (
            job_metrics[
                "cluster_id"
            ]
            .dropna()
            .astype(int)
            .value_counts()
        )

        top_cluster_ids = (
            cluster_counts
            .head(10)
            .index
            .tolist()
        )

        other_cluster_ids = [
            cluster_id
            for cluster_id in cluster_counts.index
            if cluster_id not in top_cluster_ids
        ]

        has_other = bool(
            other_cluster_ids
        )

        # Reserve category 0 for the aggregated "Other" group.
        #
        # Top individual clusters then occupy categories 1...N,
        # ordered from most to least populous.
        other_index = 0

        cluster_to_colour_index = {
            cluster_id: index
            for index, cluster_id
            in enumerate(
                top_cluster_ids,
                start=1,
            )
        }

        def cluster_colour_index(
            cluster_id
        ):
            if pd.isna(cluster_id):
                return (
                    other_index
                    if other_index is not None
                    else 0
                )

            cluster_id = int(
                cluster_id
            )

            if cluster_id in cluster_to_colour_index:
                return cluster_to_colour_index[
                    cluster_id
                ]

            return (
                other_index
                if other_index is not None
                else 0
            )

        colour_indices = (
            job_metrics[
                "cluster_id"
            ]
            .map(
                cluster_colour_index
            )
            .to_numpy(int)
        )

        # ---------------------------------------------
        # Colours
        # ---------------------------------------------

        n_top = len(
            top_cluster_ids
        )

        n_categories = (
            n_top
            + int(has_other)
        )

        base_cmap = plt.get_cmap(
            "tab10"
        )

        if has_other:

            other_index = 0

            cluster_to_colour_index = {
                cluster_id: index
                for index, cluster_id
                in enumerate(
                    top_cluster_ids,
                    start=1,
                )
            }

        else:

            cluster_to_colour_index = {
                cluster_id: index
                for index, cluster_id
                in enumerate(
                    top_cluster_ids
                )
            }

            other_index = None
        
        if has_other:
            # Other is category 0 and deliberately neutral grey.
            colours = [
                "#bdbdbd"
            ]

            colours.extend(
                base_cmap(i)
                for i in range(n_top)
            )

        else:
            colours = [
                base_cmap(i)
                for i in range(n_top)
            ]


        cmap = ListedColormap(
            colours
        )

        boundaries = np.arange(
            -0.5,
            n_categories + 0.5,
            1.0,
        )

        norm = BoundaryNorm(
            boundaries,
            cmap.N,
        )

        # ---------------------------------------------
        # Scatter
        # ---------------------------------------------

        scatter = axis.scatter(
            job_metrics[
                "Ligand RMSD"
            ],
            job_metrics[
                "Glide Density Score"
            ],
            c=colour_indices,
            cmap=cmap,
            norm=norm,
            alpha=0.85,
        )

        axis.set_title(
            job_name
        )

        axis.set_xlabel(
            "Reference RMSD (Å)"
        )

        axis.set_ylabel(
            "Glide Density Score"
        )

        # ---------------------------------------------
        # Discrete cluster legend
        # ---------------------------------------------

        legend_handles = []

        # Largest cluster first.
        for cluster_id in top_cluster_ids:

            colour_index = (
                cluster_to_colour_index[
                    cluster_id
                ]
            )

            count = int(
                cluster_counts.loc[
                    cluster_id
                ]
            )

            legend_handles.append(
                Patch(
                    facecolor=cmap(
                        norm(colour_index)
                    ),
                    edgecolor="none",
                    label=(
                        f"C{cluster_id}: "
                        f"{count} poses"
                    ),
                )
            )

        # Other always appears last.
        if has_other:

            other_pose_count = int(
                cluster_counts.loc[
                    other_cluster_ids
                ].sum()
            )

            legend_handles.append(
                Patch(
                    facecolor="#bdbdbd",
                    edgecolor="none",
                    label=(
                        f"Other: {other_pose_count} poses "
                        f"({len(other_cluster_ids)} clusters)"
                    ),
                )
            )

        axis.legend(
            handles=legend_handles,
            title="Pose clusters",
            loc="upper center",
            bbox_to_anchor=(
                0.5,
                -0.20,
            ),
            ncol=2,
            frameon=False,
            fontsize=8,
            title_fontsize=9,
        )

    # Hide unused subplot axes.
    for axis in axes.flat[n_jobs:]:
        axis.set_visible(
            False
        )

    fig.suptitle(
        "Clustered GlideEM density-score x RMSD",
        y=1.02,
    )

    fig.tight_layout(
        rect=[
            0.02,
            0.02,
            0.98,
            0.94,
        ]
    )

    fig.savefig(
        output_dir /
        "density_rmsd_small_multiples.png",
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(fig)

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
    backend.run_glide(inp, job["dir"], job["dir"] / f"{job['name']}_grid.log")
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
        "Export top-scoring MCSS poses for PyMOL visualization."
    ),
)
def export_mcss_pymol(
    backend,
    ligand,
    config,
    **kwargs,
):
    complex_job = kwargs["complex_job"]
    poseviewer = Path(kwargs["poseviewer"])
    top_n = int(kwargs.get("top_n", 5))

    cryoem_reference = Path(
        kwargs["cryoem_reference"]
    )

    glideem_reference = Path(
        complex_job["reference"]
    )

    pymol_cfg = config.get(
        "pymol",
        {}
    )
    
    binding_site_cutoff = float(
        pymol_cfg.get(
            "binding_site_cutoff",
            5.0,
        )
    )

    if not cryoem_reference.is_file():
        raise RuntimeError(
            f"Original cryo-EM reference ligand not found: "
            f"{cryoem_reference}"
        )

    if not glideem_reference.is_file():
        raise RuntimeError(
            f"GlideEM reference pose not found: "
            f"{glideem_reference}"
        )

    name = complex_job["name"]

    out_dir = (
        poseviewer.parent /
        "pymol"
    )

    receptor = None
    poses = []

    for st in structure.StructureReader(
        str(poseviewer)
    ):
        # Receptor CT in Glide PoseViewer
        if st.property.get(
            "b_glide_receptor",
            False,
        ):
            receptor = st
            continue

        score = glide_score(st)

        if score is None:
            continue

        cid = compound_id(st)

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
            f"No scored docked poses found in {poseviewer}"
        )

    out_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    poses.sort(
        key=lambda item: item[0]
    )

    top = poses[:top_n]

    # ---------------------------------------------
    # Receptor
    # ---------------------------------------------

    receptor_file = (
        out_dir /
        "receptor.pdb"
    )

    with structure.StructureWriter(
        str(receptor_file)
    ) as writer:
        writer.append(receptor)

    # ---------------------------------------------
    # Structural reference ligands
    # ---------------------------------------------

    cryoem_structure_file = (
        out_dir /
        "cryoem_structure.pdb"
    )

    cryoem_structures = list(
        structure.StructureReader(
            str(cryoem_reference)
        )
    )

    if not cryoem_structures:
        raise RuntimeError(
            f"No structure found in "
            f"{cryoem_reference}"
        )

    with structure.StructureWriter(
        str(cryoem_structure_file)
    ) as writer:
        writer.append(
            cryoem_structures[0]
        )

    glideem_reference_file = (
        out_dir /
        "glideem_reference_pose.pdb"
    )

    glideem_structures = list(
        structure.StructureReader(
            str(glideem_reference)
        )
    )

    if not glideem_structures:
        raise RuntimeError(
            f"No structure found in "
            f"{glideem_reference}"
        )

    with structure.StructureWriter(
        str(glideem_reference_file)
    ) as writer:
        writer.append(
            glideem_structures[0]
        )

    # ---------------------------------------------
    # Docked ligands
    # ---------------------------------------------

    exported = []

    for rank, (
        score,
        cid,
        st,
    ) in enumerate(
        top,
        start=1,
    ):
        # Make compound ID safe for filenames/PyMOL object names.
        safe_cid = re.sub(
            r"[^A-Za-z0-9_]+",
            "_",
            str(cid),
        ).strip("_")

        if not safe_cid:
            safe_cid = "unknown"

        object_name = (
            f"rank{rank:02d}_{safe_cid}"
        )

        pose_file = (
            out_dir /
            f"{object_name}.pdb"
        )

        with structure.StructureWriter(
            str(pose_file)
        ) as writer:
            writer.append(st)

        exported.append(
            (
                rank,
                score,
                cid,
                object_name,
                pose_file,
            )
        )

    # ---------------------------------------------
    # Write ranking table
    # ---------------------------------------------

    ranking_file = (
        out_dir /
        f"top{len(exported)}.csv"
    )

    with ranking_file.open(
        "w",
        newline="",
    ) as handle:

        writer = csv.writer(handle)

        writer.writerow(
            [
                "rank",
                "compound_id",
                "glide_score",
                "pymol_object",
                "pose_file",
            ]
        )

        for (
            rank,
            score,
            cid,
            object_name,
            pose_file,
        ) in exported:

            writer.writerow(
                [
                    rank,
                    cid,
                    score,
                    object_name,
                    pose_file.name,
                ]
            )

    # ---------------------------------------------
    # PyMOL script
    # ---------------------------------------------

    script = [
        "reinitialize",
        "bg_color white",
        "",
        "load receptor.pdb, receptor",
        "",
        "hide everything",
        "show cartoon, receptor",

        # Receptor carbons/backbone remain neutral grey.
        "color gray80, receptor",

        # Protein atoms by element.
        "color blue, (receptor and elem N)",
        "color red, (receptor and elem O)",
        "color yellow, (receptor and elem S)",
        "color orange, (receptor and elem P)",
        "color white, (receptor and elem H)",

        "",
    ]

    script.extend(
        [
            # Original experimental/cryo-EM ligand
            "load cryoem_structure.pdb, cryoem_structure",
            "show sticks, cryoem_structure",
            "color marine, "
            "(cryoem_structure and elem C)",
            "color blue, "
            "(cryoem_structure and elem N)",
            "color red, "
            "(cryoem_structure and elem O)",
            "color yellow, "
            "(cryoem_structure and elem S)",
            "color orange, "
            "(cryoem_structure and elem P)",
            "color white, "
            "(cryoem_structure and elem H)",
            "color green, "
            "(cryoem_structure and elem F)",
            "color green, "
            "(cryoem_structure and elem Cl)",
            "color brown, "
            "(cryoem_structure and elem Br)",
            "color purple, "
            "(cryoem_structure and elem I)",
            "",

            # GlideEM hypothesis used as the MCSS reference
            "load glideem_reference_pose.pdb, "
            "glideem_reference_pose",
            "show sticks, glideem_reference_pose",
            "color salmon, "
            "(glideem_reference_pose and elem C)",
            "color blue, "
            "(glideem_reference_pose and elem N)",
            "color red, "
            "(glideem_reference_pose and elem O)",
            "color yellow, "
            "(glideem_reference_pose and elem S)",
            "color orange, "
            "(glideem_reference_pose and elem P)",
            "color white, "
            "(glideem_reference_pose and elem H)",
            "color green, "
            "(glideem_reference_pose and elem F)",
            "color green, "
            "(glideem_reference_pose and elem Cl)",
            "color brown, "
            "(glideem_reference_pose and elem Br)",
            "color purple, "
            "(glideem_reference_pose and elem I)",
            "",
        ]
    )

    # Each ligand gets a different carbon colour.
    # Heteroatoms are subsequently recoloured by element.
    carbon_colours = [
        "cyan",
        "green",
        "yellow",
        "orange",
        "magenta",
    ]

    ligand_objects = []

    for (
        rank,
        score,
        cid,
        object_name,
        pose_file,
    ) in exported:

        carbon_colour = carbon_colours[
            (rank - 1) % len(carbon_colours)
        ]

        ligand_objects.append(
            object_name
        )

        script.extend(
            [
                f"load {pose_file.name}, {object_name}",
                f"show sticks, {object_name}",

                # Carbon distinguishes individual ligands.
                f"color {carbon_colour}, "
                f"({object_name} and elem C)",

                # Conventional element colours.
                f"color blue, "
                f"({object_name} and elem N)",

                f"color red, "
                f"({object_name} and elem O)",

                f"color yellow, "
                f"({object_name} and elem S)",

                f"color orange, "
                f"({object_name} and elem P)",

                f"color white, "
                f"({object_name} and elem H)",

                # Halogens
                f"color green, "
                f"({object_name} and elem F)",

                f"color green, "
                f"({object_name} and elem Cl)",

                f"color brown, "
                f"({object_name} and elem Br)",

                f"color purple, "
                f"({object_name} and elem I)",

                "",
            ]
        )

    # Group the poses in the PyMOL object tree.
    if ligand_objects:
        script.append(
            "group top_mcss_poses, "
            + " ".join(ligand_objects)
        )

    # ---------------------------------------------
    # Binding-site protein residues
    # ---------------------------------------------

    ligand_selection = (
        "cryoem_structure or "
        "glideem_reference_pose or "
        "top_mcss_poses"
    )

    script.extend(
        [
            "",
            "# Protein residues surrounding the displayed ligands",
            (
                "select binding_site, "
                f"byres (receptor within "
                f"{binding_site_cutoff:.1f} of "
                f"({ligand_selection}))"
            ),
            "show sticks, binding_site",

            # Protein carbon remains neutral grey.
            "color gray80, "
            "(binding_site and elem C)",

            # Standard heteroatom colouring.
            "color blue, "
            "(binding_site and elem N)",

            "color red, "
            "(binding_site and elem O)",

            "color yellow, "
            "(binding_site and elem S)",

            "color orange, "
            "(binding_site and elem P)",

            "color white, "
            "(binding_site and elem H)",

            "",
        ]
    )

    script.extend(
        [
            "",
            "orient top_mcss_poses",
            "zoom top_mcss_poses, 8",
            "",
        ]
    )

    script_file = (
        out_dir /
        f"pymol_all_mcss_poses_cryo_glideem_ref.pml"
    )

    script_file.write_text(
        "\n".join(script)
        + "\n"
    )

    logger.debug(
        "[%s] Exported top %d MCSS poses "
        "for PyMOL to %s",
        name,
        len(exported),
        out_dir,
    )

    return script_file