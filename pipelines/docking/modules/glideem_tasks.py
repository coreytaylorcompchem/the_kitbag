from __future__ import annotations

import shutil
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

from modules.utils.mcss_docking_helpers import (
    export_mcss_pose_bundle,
    safe_pymol_name
)

from modules.utils.prime_mmgbsa_helpers import (
    _write_mmgbsa_input,
    _parse_mmgbsa_output,
)

from modules.utils.glideem_campaign_analysis import (
    cluster_global_pose_families,
    centroid_dispersion_summary,
    density_difference,
    density_variability,
    family_medoid_indices,
    family_support_summary,
    load_aligned_glideem_poses,
    make_common_grid,
    pooled_probability_density,
    pose_probability_density,
    read_first_structure,
    read_reference_ligand_coordinates,
    save_global_family_dendrogram,
    save_grid_metadata,
    save_numpy_density,
    save_structure_family_heatmap,
    save_structure_js_heatmap,
    structure_ambiguity_summary,
    structure_family_composition,
    structure_js_distance_matrix,
    write_aligned_pose_archive,
    transform_structure_copy,
    extract_ca_chains,
    select_best_sequence_matched_chain_pair,
    sequence_mapped_pocket_coordinates,
    robust_kabsch_transform,
    pairwise_centroid_distance_matrix,
    build_mobile_to_reference_position_map,
    receptor_contact_atoms,
    assign_second_level_subfamilies,
    pose_contact_fingerprint,
    ligand_orientation_descriptors,
    export_second_level_pymol,
    save_second_level_family_heatmap,
    second_level_contact_summary,
)

from modules.utils.glideem_ensemble_support import (
    run_ensemble_support_analysis,
)

from modules.utils.glideem_pocket_evidence import (
    run_pocket_evidence_analysis,
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

    output_root = (
        Path(
            config[
                "output_dir"
            ]
        )
        .expanduser()
        .resolve()
    )

    output_dir = (
        output_root
        / "ligprep"
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    ligprep_cfg = config.get(
        "ligprep",
        {},
    )

    already_prepared = bool(
        ligprep_cfg.get(
            "already_prepared",
            False,
        )
    )

    output_file = (
        output_dir
        / "activity_ligands_ligprep.maegz"
    )

    source_record = (
        output_dir
        / "activity_ligands_source.txt"
    )

    # --------------------------------------------------
    # Already-prepared ligand library
    #
    # Copy it into the canonical pipeline location without
    # running LigPrep. Downstream tasks therefore continue
    # to consume exactly the same filename regardless of
    # preparation route.
    # --------------------------------------------------

    if already_prepared:
        allowed_prepared_suffixes = {
            ".mae",
            ".maegz",
        }

        if (
            input_suffix
            not in allowed_prepared_suffixes
        ):
            raise ValueError(
                "ligprep.already_prepared=true requires "
                "a prepared Maestro ligand library "
                "(.mae or .maegz); received: "
                f"{input_path}"
            )

        can_reuse = (
            output_file.is_file()
            and output_file.stat().st_size > 0
            and source_record.is_file()
            and source_record.read_text().strip()
            == str(
                input_path
            )
            and not config.get(
                "force",
                False,
            )
        )

        if can_reuse:
            logger.info(
                "Using existing copied prepared ligand "
                "library: %s",
                output_file,
            )

            return output_file

        logger.info(
            "Using pre-prepared activity ligand library; "
            "LigPrep will not be run: %s",
            input_path,
        )

        if (
            input_path
            == output_file
        ):
            logger.info(
                "Prepared ligand library is already at "
                "the canonical pipeline location: %s",
                output_file,
            )

        else:
            shutil.copy2(
                input_path,
                output_file,
            )

            logger.info(
                "Copied prepared ligand library: "
                "%s -> %s",
                input_path,
                output_file,
            )

        if (
            not output_file.is_file()
            or output_file.stat().st_size == 0
        ):
            raise RuntimeError(
                "Failed to copy prepared activity ligand "
                f"library to: {output_file}"
            )

        source_record.write_text(
            str(
                input_path
            )
            + "\n"
        )

        return output_file

    # --------------------------------------------------
    # Normal LigPrep route
    # --------------------------------------------------

    allowed_unprepared_suffixes = {
        ".sdf",
        ".sd",
        ".csv",
    }

    if (
        input_suffix
        not in allowed_unprepared_suffixes
    ):
        raise ValueError(
            "activity_ligands must be an SDF/SD or CSV "
            "when LigPrep preparation is enabled; "
            f"received: {input_path}"
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

    # output_file = (
    #     output_dir /
    #     "activity_ligands_ligprep.maegz"
    # )

    # # Reuse existing output unless explicitly forcing rerun

    # source_record = (
    #     output_dir
    #     / "activity_ligands_source.txt"
    # )

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

    # ligprep_cfg = config.get(
    #     "ligprep",
    #     {}
    # )

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

    source_record.write_text(
        str(
            input_path
        )
        + "\n"
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
        # Required by campaign-level landscape analysis
        "archive": archive,
        "prepared_receptor": prepared_receptor,
        "cryoem_reference": cryoem_reference,
        "reference_ligand_resname": ligand_resname,
        "ligand_from_smiles": ligand_from_smiles,
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
    "calculate_mmgbsa",
    category="Prime",
    description=(
        "Rescore MCSS-constrained Glide poses "
        "using Prime MM-GBSA."
    ),
)
def calculate_mmgbsa(
    backend,
    ligand,
    config,
    **kwargs,
):
    complex_job = kwargs[
        "complex_job"
    ]

    poseviewer = (
        Path(kwargs["poseviewer"])
        .expanduser()
        .resolve()
    )

    if not poseviewer.is_file():
        raise FileNotFoundError(
            f"MCSS PoseViewer not found: "
            f"{poseviewer}"
        )

    if poseviewer.stat().st_size == 0:
        raise RuntimeError(
            f"MCSS PoseViewer is empty: "
            f"{poseviewer}"
        )

    mmgbsa_cfg = config.get(
        "mmgbsa",
        {},
    )

    job_name = (
        f"{complex_job['name']}"
        "_mmgbsa"
    )

    output_dir = (
        Path(kwargs["mmgbsa_root"])
        .expanduser()
        .resolve()
        / complex_job["name"]
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    input_file = (
        output_dir
        / f"{job_name}_input.maegz"
    )

    input_manifest_file = (
        output_dir
        / f"{job_name}_input_poses.csv"
    )

    results_csv = (
        output_dir
        / f"{job_name}_results.csv"
    )

    output_file = (
        output_dir
        / f"{job_name}-out.maegz"
    )

    launcher_log = (
        output_dir
        / f"{job_name}_launcher.log"
    )

    can_reuse = (
        results_csv.is_file()
        and results_csv.stat().st_size > 0
        and not config.get(
            "force",
            False,
        )
    )

    if can_reuse:
        logger.info(
            "[%s] Using existing MM-GBSA "
            "results: %s",
            complex_job["name"],
            results_csv,
        )

        return results_csv

    input_manifest = (
        _write_mmgbsa_input(
            poseviewer=poseviewer,
            output_file=input_file,
        )
    )

    input_manifest.to_csv(
        input_manifest_file,
        index=False,
    )

    logger.info(
        "[%s] Running Prime MM-GBSA for "
        "%d MCSS poses",
        complex_job["name"],
        len(input_manifest),
    )

    backend.run_mmgbsa(
        input_file=input_file,
        cwd=output_dir,
        log_path=launcher_log,
        job_name=job_name,
        host=mmgbsa_cfg.get(
            "host",
            backend.host,
        ),
    )

    if not output_file.is_file():
        candidates = sorted(
            output_dir.glob(
                f"{job_name}*out*.maegz"
            )
        )

        if len(candidates) == 1:
            output_file = candidates[
                0
            ]

        else:
            raise RuntimeError(
                "Prime MM-GBSA completed but the "
                "expected output structure file was "
                f"not found: {output_file}. "
                f"Candidates: {candidates}"
            )

    if output_file.stat().st_size == 0:
        raise RuntimeError(
            "Prime MM-GBSA output is empty: "
            f"{output_file}"
        )

    results = _parse_mmgbsa_output(
        output_file=output_file,
        input_manifest=input_manifest,
    )

    results.insert(
        0,
        "dock_name",
        complex_job["name"],
    )

    results.insert(
        1,
        "source_name",
        complex_job[
            "source_name"
        ],
    )

    results.to_csv(
        results_csv,
        index=False,
    )

    logger.info(
        "[%s] Prime MM-GBSA completed: "
        "%d/%d poses have dG Bind values",
        complex_job["name"],
        int(
            results[
                "mmgbsa_dg_bind"
            ].notna().sum()
        ),
        len(results),
    )

    return results_csv

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

    output_root = (
        Path(
            kwargs[
                "output_dir"
            ]
        )
        .expanduser()
        .resolve()
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

    name = str(
        complex_job[
            "name"
        ]
    )

    source_name = str(
        complex_job[
            "source_name"
        ]
    )

    # complex_job["name"] is expected to look like:
    #
    #     structure_51052_pose01
    #
    # while source_name is:
    #
    #     structure_51052
    #
    # Remove the source prefix to obtain:
    #
    #     pose01
    prefix = (
        source_name
        + "_"
    )

    if name.startswith(
        prefix
    ):
        pose_name = name[
            len(
                prefix
            ):
        ]

    else:
        # Defensive fallback in case a future complex-job naming
        # convention does not begin with the source structure name.
        pose_name = name

    pose_name = safe_pymol_name(
        pose_name,
        fallback="pose",
    )

    out_dir = (
        output_root
        / source_name
        / f"{pose_name}_pymol"
    )

    out_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    logger.debug(
        "[%s] MCSS PyMOL bundle directory: %s",
        name,
        out_dir,
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

@register_task(
    "analyse_glideem_campaign_landscape",
    category="GlideEM",
    description=(
        "Align all GlideEM ensembles into a common receptor frame "
        "and analyse the campaign-level 3D pose landscape."
    ),
)
def analyse_glideem_campaign_landscape(
    backend,
    ligand,
    config,
    **kwargs,
):
    analysis_results = list(
        kwargs["analysis_results"]
    )

    output_dir = (
        Path(
            kwargs["output_dir"]
        )
        .expanduser()
        .resolve()
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not analysis_results:
        raise RuntimeError(
            "No per-structure GlideEM analyses were supplied "
            "to campaign landscape analysis."
        )

    campaign_cfg = config.get(
        "glideem_campaign_analysis",
        {},
    )

    alignment_cfg = campaign_cfg.get(
        "alignment",
        {},
    )

    density_cfg = campaign_cfg.get(
        "density",
        {},
    )

    clustering_cfg = campaign_cfg.get(
        "clustering",
        {},
    )

    export_cfg = campaign_cfg.get(
        "export",
        {},
    )

    ensemble_support_cfg = campaign_cfg.get(
        "ensemble_support",
        {},
    )

    ensemble_support_enabled = bool(
        ensemble_support_cfg.get(
            "enabled",
            True,
        )
    )

    second_level_cfg = campaign_cfg.get(
        "second_level",
        {},
    )    

    pocket_cutoff = float(
        alignment_cfg.get(
            "pocket_cutoff",
            10.0,
        )
    )

    alignment_atom_name = str(
        alignment_cfg.get(
            "atom_name",
            "CA",
        )
    ).strip().upper()

    minimum_common_atoms = int(
        alignment_cfg.get(
            "minimum_common_atoms",
            20,
        )
    )

    minimum_chain_length = int(
        alignment_cfg.get(
            "minimum_chain_length",
            100,
        )
    )

    minimum_sequence_identity = float(
        alignment_cfg.get(
            "minimum_sequence_identity",
            0.70,
        )
    )

    minimum_reference_coverage = float(
        alignment_cfg.get(
            "minimum_reference_coverage",
            0.70,
        )
    )

    coordinate_outlier_cutoff = float(
        alignment_cfg.get(
            "coordinate_outlier_cutoff",
            3.0,
        )
    )

    requested_reference_name = (
        alignment_cfg.get(
            "reference_structure"
        )
    )

    spacing = float(
        density_cfg.get(
            "spacing",
            0.75,
        )
    )

    padding = float(
        density_cfg.get(
            "padding",
            4.0,
        )
    )

    gaussian_sigma = float(
        density_cfg.get(
            "gaussian_sigma",
            2.0,
        )
    )

    # distance_metric = str(
    #     clustering_cfg.get(
    #         "distance_metric",
    #         "jensenshannon",
    #     )
    # ).strip().lower()

    linkage_method = str(
        clustering_cfg.get(
            "linkage_method",
            "average",
        )
    ).strip().lower()

    centroid_distance_cutoff = float(
        clustering_cfg.get(
            "centroid_distance_cutoff",
            2.5,
        )
    )

    if centroid_distance_cutoff <= 0.0:
        raise ValueError(
            "glideem_campaign_analysis.clustering."
            "centroid_distance_cutoff must be positive."
        )

    contact_cutoff = float(
        second_level_cfg.get(
            "contact_cutoff",
            4.5,
        )
    )

    polar_contact_cutoff = float(
        second_level_cfg.get(
            "polar_contact_cutoff",
            3.5,
        )
    )

    orientation_weight = float(
        second_level_cfg.get(
            "orientation_weight",
            0.40,
        )
    )

    contact_weight = float(
        second_level_cfg.get(
            "contact_weight",
            0.60,
        )
    )

    subfamily_distance_cutoff = float(
        second_level_cfg.get(
            "subfamily_distance_cutoff",
            0.45,
        )
    )

    minimum_axis_anisotropy = float(
        second_level_cfg.get(
            "minimum_axis_anisotropy",
            0.15,
        )
    )

    export_second_level_pymol_enabled = bool(
        second_level_cfg.get(
            "export_pymol",
            True,
        )
    )

    if contact_cutoff <= 0.0:
        raise ValueError(
            "second_level.contact_cutoff must be positive."
        )

    if (
        polar_contact_cutoff <= 0.0
        or polar_contact_cutoff > contact_cutoff
    ):
        raise ValueError(
            "second_level.polar_contact_cutoff must be positive "
            "and no greater than contact_cutoff."
        )

    if (
        orientation_weight < 0.0
        or contact_weight < 0.0
        or (
            orientation_weight
            + contact_weight
        ) <= 0.0
    ):
        raise ValueError(
            "Second-level weights must be non-negative and "
            "must have a positive sum."
        )

    if spacing <= 0.0:
        raise ValueError(
            "glideem_campaign_analysis.density.spacing "
            "must be positive."
        )

    if gaussian_sigma < 0.0:
        raise ValueError(
            "glideem_campaign_analysis.density.gaussian_sigma "
            "cannot be negative."
        )

    pocket_evidence_cfg = campaign_cfg.get(
        "pocket_evidence",
        {},
    )

    pocket_evidence_enabled = bool(
        pocket_evidence_cfg.get(
            "enabled",
            True,
        )
    )

    # --------------------------------------------------
    # Validate and index source analyses
    # --------------------------------------------------

    valid_results = []

    for result in analysis_results:
        required_fields = {
            "job_name",
            "archive",
            "prepared_receptor",
            "cryoem_reference",
            "reference_ligand_resname",
            "ligand_from_smiles",
            "metrics",
        }

        missing_fields = (
            required_fields
            - set(
                result
            )
        )

        if missing_fields:
            logger.warning(
                "[%s] Excluded from campaign landscape because "
                "analysis result fields are missing: %s",
                result.get(
                    "job_name",
                    "unknown",
                ),
                sorted(
                    missing_fields
                ),
            )
            continue

        archive = Path(
            result[
                "archive"
            ]
        )

        receptor_file = Path(
            result[
                "prepared_receptor"
            ]
        )

        reference_file = Path(
            result[
                "cryoem_reference"
            ]
        )

        metrics_file = Path(
            result[
                "metrics"
            ]
        )

        missing_files = [
            path
            for path in (
                archive,
                receptor_file,
                reference_file,
                metrics_file,
            )
            if not path.is_file()
        ]

        if missing_files:
            logger.warning(
                "[%s] Excluded from campaign landscape because "
                "files are missing: %s",
                result[
                    "job_name"
                ],
                missing_files,
            )
            continue

        valid_results.append(
            result
        )

    if len(
        valid_results
    ) < 2:
        raise RuntimeError(
            "Campaign landscape analysis requires at least "
            "two valid GlideEM structures."
        )

    # --------------------------------------------------
    # Choose common coordinate frame
    # --------------------------------------------------

    if requested_reference_name:
        matching = [
            result
            for result in valid_results
            if result[
                "job_name"
            ] == requested_reference_name
        ]

        if not matching:
            raise ValueError(
                "Requested campaign alignment reference was "
                f"not found: {requested_reference_name}"
            )

        reference_result = matching[
            0
        ]

    else:
        reference_result = valid_results[
            0
        ]

    reference_name = reference_result[
        "job_name"
    ]

    reference_receptor = (
        read_first_structure(
            Path(
                reference_result[
                    "prepared_receptor"
                ]
            )
        )
    )

    reference_ligand_coordinates = (
        read_reference_ligand_coordinates(
            Path(
                reference_result[
                    "cryoem_reference"
                ]
            )
        )
    )

    # reference_pocket_atoms = (
    #     pocket_alignment_coordinates(
    #         receptor_st=reference_receptor,
    #         reference_ligand_coordinates=(
    #             reference_ligand_coordinates
    #         ),
    #         pocket_cutoff=pocket_cutoff,
    #         atom_name=alignment_atom_name,
    #     )
    # )

    # if (
    #     len(
    #         reference_pocket_atoms
    #     )
    #     < minimum_common_atoms
    # ):
    #     raise RuntimeError(
    #         f"Reference structure {reference_name} contains only "
    #         f"{len(reference_pocket_atoms)} pocket alignment atoms; "
    #         f"at least {minimum_common_atoms} are required."
    #     )

    reference_chains = extract_ca_chains(
        receptor_st=reference_receptor,
        atom_name=alignment_atom_name,
        minimum_chain_length=minimum_chain_length,
    )

    if not reference_chains:
        raise RuntimeError(
            f"No plausible receptor chains were found in "
            f"reference structure {reference_name}."
        )

    logger.info(
        "Campaign alignment reference: %s with %d "
        "candidate protein chains",
        reference_name,
        len(
            reference_chains
        ),
    )
    # --------------------------------------------------
    # Align every structure and load all poses
    # --------------------------------------------------

    all_pose_records = []
    alignment_rows = []
    accepted_results = []

    # Aligned receptor structures used for visual QC of the
    # campaign coordinate-frame transformation.
    aligned_receptor_records = []

    contact_receptor_atoms_by_source = {}
    campaign_reference_chain_id = None

    for result in valid_results:
        source_name = result[
            "job_name"
        ]

        receptor_st = read_first_structure(
            Path(
                result[
                    "prepared_receptor"
                ]
            )
        )

        # source_reference_coordinates = (
        #     read_reference_ligand_coordinates(
        #         Path(
        #             result[
        #                 "cryoem_reference"
        #             ]
        #         )
        #     )
        # )

        # source_pocket_atoms = (
        #     pocket_alignment_coordinates(
        #         receptor_st=receptor_st,
        #         reference_ligand_coordinates=(
        #             source_reference_coordinates
        #         ),
        #         pocket_cutoff=pocket_cutoff,
        #         atom_name=alignment_atom_name,
        #     )
        # )

        # common_keys = sorted(
        #     set(
        #         reference_pocket_atoms
        #     )
        #     & set(
        #         source_pocket_atoms
        #     )
        # )

        # if (
        #     len(
        #         common_keys
        #     )
        #     < minimum_common_atoms
        # ):
        #     logger.warning(
        #         "[%s] Excluded from campaign analysis: only %d "
        #         "common pocket alignment atoms were found",
        #         source_name,
        #         len(
        #             common_keys
        #         ),
        #     )
        #     continue

        # mobile_coordinates = np.asarray(
        #     [
        #         source_pocket_atoms[
        #             key
        #         ]
        #         for key in common_keys
        #     ],
        #     dtype=float,
        # )

        # target_coordinates = np.asarray(
        #     [
        #         reference_pocket_atoms[
        #             key
        #         ]
        #         for key in common_keys
        #     ],
        #     dtype=float,
        # )

        # rotation, translation = (
        #     kabsch_transform(
        #         mobile_coordinates,
        #         target_coordinates,
        #     )
        # )

        # fit_rmsd = alignment_rmsd(
        #     mobile_coordinates=mobile_coordinates,
        #     reference_coordinates=target_coordinates,
        #     rotation=rotation,
        #     translation=translation,
        # )

        mobile_chains = extract_ca_chains(
            receptor_st=receptor_st,
            atom_name=alignment_atom_name,
            minimum_chain_length=minimum_chain_length,
        )

        if not mobile_chains:
            logger.warning(
                "[%s] Excluded from campaign analysis: "
                "no plausible protein chains were found",
                source_name,
            )
            continue

        try:
            chain_match = (
                select_best_sequence_matched_chain_pair(
                    reference_chains=reference_chains,
                    mobile_chains=mobile_chains,
                    minimum_identity=(
                        minimum_sequence_identity
                    ),
                    minimum_reference_coverage=(
                        minimum_reference_coverage
                    ),
                )
            )

            pocket_mapping = (
                sequence_mapped_pocket_coordinates(
                    chain_match=chain_match,
                    reference_ligand_coordinates=(
                        reference_ligand_coordinates
                    ),
                    pocket_cutoff=pocket_cutoff,
                )
            )

            robust_fit = robust_kabsch_transform(
                mobile_coordinates=(
                    pocket_mapping[
                        "mobile_coordinates"
                    ]
                ),
                reference_coordinates=(
                    pocket_mapping[
                        "reference_coordinates"
                    ]
                ),
                minimum_atoms=minimum_common_atoms,
                outlier_cutoff=(
                    coordinate_outlier_cutoff
                ),
            )

        except Exception as exc:
            logger.warning(
                "[%s] Excluded from campaign alignment: %s",
                source_name,
                exc,
            )
            continue

        rotation = robust_fit[
            "rotation"
        ]

        translation = robust_fit[
            "translation"
        ]

        fit_rmsd = robust_fit[
            "fit_rmsd"
        ]

        if source_name == reference_name:
            campaign_reference_chain_id = (
                chain_match[
                    "reference_chain_id"
                ]
            )

        retained_mask = robust_fit[
            "retained_mask"
        ]

        retained_residue_pairs = [
            residue_pair
            for residue_pair, retained in zip(
                pocket_mapping[
                    "residue_pairs"
                ],
                retained_mask,
            )
            if retained
        ]

        # ----------------------------------------------
        # Apply exactly the same rigid transformation to
        # the complete receptor for alignment QC.
        # ----------------------------------------------

        receptor_for_export = receptor_st

        if not bool(
            result["ligand_from_smiles"]
        ):
            ligand_indices = find_glideem_ligand_indices(
                receptor_st,
                ligand_resname=(
                    result[
                        "reference_ligand_resname"
                    ]
                ),
                ligand_from_smiles=False,
            )

            if ligand_indices:
                ligand_index_set = set(
                    ligand_indices
                )

                receptor_atom_indices = [
                    atom.index
                    for atom in receptor_st.atom
                    if atom.index not in ligand_index_set
                ]

                receptor_for_export = (
                    receptor_st.extract(
                        receptor_atom_indices
                    )
                )

        aligned_receptor = transform_structure_copy(
            input_structure=receptor_for_export,
            rotation=rotation,
            translation=translation,
        )


        aligned_receptor.title = (
            f"{source_name}_aligned_to_{reference_name}"
        )

        aligned_receptor.property[
            "s_user_source_name"
        ] = source_name

        aligned_receptor.property[
            "s_user_alignment_reference"
        ] = reference_name

        aligned_receptor.property[
            "i_user_common_alignment_atoms"
        ] = int(
            robust_fit["n_retained_atoms"]
        )

        aligned_receptor.property[
            "r_user_pocket_alignment_rmsd"
        ] = float(
            fit_rmsd
        )

        aligned_receptor_records.append(
            {
                "source_name": source_name,
                "structure": aligned_receptor,
            }
        )

        mobile_to_reference_position = (
            build_mobile_to_reference_position_map(
                chain_match
            )
        )

        contact_receptor_atoms_by_source[
            source_name
        ] = receptor_contact_atoms(
            aligned_receptor=aligned_receptor,
            mobile_to_reference_position=(
                mobile_to_reference_position
            ),
        )

        if not contact_receptor_atoms_by_source[
            source_name
        ]:
            raise RuntimeError(
                f"[{source_name}] No sequence-mapped receptor "
                "contact atoms were generated."
            )        

        aligned_poses = (
            load_aligned_glideem_poses(
                archive=Path(
                    result[
                        "archive"
                    ]
                ),
                ligand_resname=result[
                    "reference_ligand_resname"
                ],
                ligand_from_smiles=bool(
                    result[
                        "ligand_from_smiles"
                    ]
                ),
                rotation=rotation,
                translation=translation,
            )
        )

        if not aligned_poses:
            logger.warning(
                "[%s] No usable aligned ligand poses; excluded",
                source_name,
            )
            continue

        metrics = pd.read_csv(
            result[
                "metrics"
            ]
        )

        metrics_by_pose = (
            metrics
            .set_index(
                "pose_number"
            )
            .to_dict(
                orient="index"
            )
        )

        accepted_results.append(
            result
        )

        alignment_rows.append(
            {
                "source_name": source_name,
                "alignment_reference": reference_name,
                "reference_chain": (
                    chain_match[
                        "reference_chain_id"
                    ]
                ),
                "mobile_chain": (
                    chain_match[
                        "mobile_chain_id"
                    ]
                ),
                "reference_chain_length": len(
                    chain_match[
                        "reference_residues"
                    ]
                ),
                "mobile_chain_length": len(
                    chain_match[
                        "mobile_residues"
                    ]
                ),
                "sequence_identity": float(
                    chain_match[
                        "alignment"
                    ][
                        "identity"
                    ]
                ),
                "reference_sequence_coverage": float(
                    chain_match[
                        "alignment"
                    ][
                        "first_coverage"
                    ]
                ),
                "mobile_sequence_coverage": float(
                    chain_match[
                        "alignment"
                    ][
                        "second_coverage"
                    ]
                ),
                "n_sequence_aligned_residues": int(
                    chain_match[
                        "alignment"
                    ][
                        "aligned_pair_count"
                    ]
                ),
                "n_initial_pocket_atoms": int(
                    robust_fit[
                        "n_initial_atoms"
                    ]
                ),
                "n_retained_pocket_atoms": int(
                    robust_fit[
                        "n_retained_atoms"
                    ]
                ),
                "n_rejected_pocket_atoms": int(
                    robust_fit[
                        "n_initial_atoms"
                    ]
                    - robust_fit[
                        "n_retained_atoms"
                    ]
                ),
                "pocket_alignment_rmsd": float(
                    fit_rmsd
                ),
                "median_final_pocket_residual": float(
                    np.median(
                        robust_fit[
                            "final_all_residuals"
                        ][
                            retained_mask
                        ]
                    )
                ),
                "max_retained_pocket_residual": float(
                    np.max(
                        robust_fit[
                            "final_all_residuals"
                        ][
                            retained_mask
                        ]
                    )
                ),
                "n_aligned_poses": len(
                    aligned_poses
                ),
            }
        )

        for pose in aligned_poses:
            pose_number = int(
                pose[
                    "pose_number"
                ]
            )

            metric_record = metrics_by_pose.get(
                pose_number,
                {},
            )

            all_pose_records.append(
                {
                    "source_name": source_name,
                    "pose_number": pose_number,
                    "coordinates": pose[
                        "coordinates"
                    ],
                    "elements": pose[
                        "elements"
                    ],
                    "structure": pose[
                        "structure"
                    ],
                    "centroid": pose[
                        "centroid"
                    ],
                    "local_cluster_id": metric_record.get(
                        "cluster_id",
                        np.nan,
                    ),
                    "reference_rmsd": metric_record.get(
                        "Ligand RMSD",
                        np.nan,
                    ),
                    "density_score": metric_record.get(
                        "Glide Density Score",
                        np.nan,
                    ),
                    "docking_score": metric_record.get(
                        "Glide Docking Score",
                        np.nan,
                    ),
                }
            )

        logger.info(
            "[%s] Aligned %d poses using sequence-matched "
            "chains %r -> %r; identity=%.3f, coverage=%.3f, "
            "pocket atoms=%d/%d, RMSD=%.3f A",
            source_name,
            len(
                aligned_poses
            ),
            chain_match[
                "reference_chain_id"
            ],
            chain_match[
                "mobile_chain_id"
            ],
            chain_match[
                "alignment"
            ][
                "identity"
            ],
            chain_match[
                "alignment"
            ][
                "first_coverage"
            ],
            robust_fit[
                "n_retained_atoms"
            ],
            robust_fit[
                "n_initial_atoms"
            ],
            fit_rmsd,
        )

    if len(
        accepted_results
    ) < 2:
        raise RuntimeError(
            "Fewer than two structures passed campaign alignment."
        )

    if len(
        all_pose_records
    ) < 2:
        raise RuntimeError(
            "Fewer than two aligned poses are available."
        )

    if campaign_reference_chain_id is None:
        raise RuntimeError(
            "The selected reference receptor chain was not recorded "
            f"for campaign reference '{reference_name}'."
        )    

    alignment_df = pd.DataFrame(
        alignment_rows
    )

    alignment_df.to_csv(
        output_dir
        / "receptor_pocket_alignment.csv",
        index=False,
    )

    # --------------------------------------------------
    # Export aligned receptors for visual alignment QC
    # --------------------------------------------------

    aligned_receptor_file = (
        output_dir
        / "all_aligned_receptors.maegz"
    )

    with structure.StructureWriter(
        str(
            aligned_receptor_file
        )
    ) as writer:
        for record in aligned_receptor_records:
            writer.append(
                record[
                    "structure"
                ]
            )

    logger.info(
        "Exported %d aligned receptors for coordinate-frame "
        "QC: %s",
        len(
            aligned_receptor_records
        ),
        aligned_receptor_file,
    )

    # --------------------------------------------------
    # Build one common grid
    # --------------------------------------------------

    common_grid = make_common_grid(
        all_coordinates=[
            record[
                "coordinates"
            ]
            for record in all_pose_records
        ],
        spacing=spacing,
        padding=padding,
    )

    save_grid_metadata(
        common_grid,
        output_dir
        / "campaign_density_grid.json",
    )

    logger.info(
        "Campaign density grid: shape=%s, spacing=%.3f A",
        common_grid[
            "shape"
        ],
        common_grid[
            "spacing"
        ],
    )

        # --------------------------------------------------
    # One normalized 3D density per pose
    #
    # These densities are retained for:
    #   - structure-level spatial probability fields
    #   - family-level spatial probability fields
    #   - structure-to-structure JS comparisons
    #
    # They are no longer used to define the first-level
    # global pose families.
    # --------------------------------------------------

    pose_densities = []

    for record in all_pose_records:
        density = pose_probability_density(
            coordinates=record[
                "coordinates"
            ],
            grid=common_grid,
            gaussian_sigma_angstrom=(
                gaussian_sigma
            ),
        )

        record[
            "density"
        ] = density

        pose_densities.append(
            density
        )

    # --------------------------------------------------
    # First-level global spatial pose families
    #
    # These families describe which broad region of the
    # aligned orthosteric pocket each ligand pose occupies.
    #
    # Distances are Euclidean distances between aligned
    # ligand centroids, expressed in Angstrom.
    # --------------------------------------------------

    global_distance_matrix = (
        pairwise_centroid_distance_matrix(
            pose_records=all_pose_records
        )
    )

    family_ids, linkage_matrix = (
        cluster_global_pose_families(
            distance_matrix=(
                global_distance_matrix
            ),
            distance_cutoff=(
                centroid_distance_cutoff
            ),
            linkage_method=linkage_method,
        )
    )

    # Attach global spatial-family assignments to every pose.
    for record, family_id in zip(
        all_pose_records,
        family_ids,
    ):
        record[
            "global_family_id"
        ] = int(
            family_id
        )

        # Explicit alias for the scientific meaning of this level.
        record[
            "global_spatial_family_id"
        ] = int(
            family_id
        )

    # Select the pose whose centroid has the lowest average distance
    # to the other members of each spatial family.
    medoid_indices = family_medoid_indices(
        distance_matrix=(
            global_distance_matrix
        ),
        family_ids=family_ids,
    )

    for record in all_pose_records:
        record[
            "is_family_medoid"
        ] = False

    for family_id, medoid_index in (
        medoid_indices.items()
    ):
        all_pose_records[
            medoid_index
        ][
            "is_family_medoid"
        ] = True

    # Persist the pairwise centroid-distance matrix.
    pose_labels = [
        (
            f"{record['source_name']}:"
            f"pose{record['pose_number']:03d}"
        )
        for record in all_pose_records
    ]

    pd.DataFrame(
        global_distance_matrix,
        index=pose_labels,
        columns=pose_labels,
    ).to_csv(
        output_dir
        / "global_pose_centroid_distances.csv"
    )

    # Persist the corresponding hierarchical clustering.
    save_global_family_dendrogram(
        linkage_matrix=linkage_matrix,
        labels=pose_labels,
        cutoff=centroid_distance_cutoff,
        output_file=(
            output_dir
            / "global_spatial_family_dendrogram.png"
        ),
        distance_label=(
            "Aligned ligand-centroid distance (A)"
        ),
        title=(
            "Global GlideEM spatial-family clustering"
        ),
    )

    logger.info(
        "Assigned %d aligned poses to %d global spatial families "
        "using a centroid-distance cutoff of %.2f A",
        len(
            all_pose_records
        ),
        len(
            set(
                family_ids.tolist()
            )
        ),
        centroid_distance_cutoff,
    )

    # --------------------------------------------------
    # Second-level orientation/contact descriptors
    # --------------------------------------------------

    for record in all_pose_records:
        orientation = ligand_orientation_descriptors(
            record[
                "coordinates"
            ]
        )

        record[
            "principal_axis"
        ] = orientation[
            "principal_axis"
        ]

        record[
            "axis_anisotropy"
        ] = orientation[
            "axis_anisotropy"
        ]

        record[
            "radius_of_gyration"
        ] = orientation[
            "radius_of_gyration"
        ]

        receptor_atoms = (
            contact_receptor_atoms_by_source[
                record[
                    "source_name"
                ]
            ]
        )

        record[
            "contact_fingerprint"
        ] = pose_contact_fingerprint(
            ligand_coordinates=record[
                "coordinates"
            ],
            ligand_elements=record[
                "elements"
            ],
            receptor_atoms=receptor_atoms,
            contact_cutoff=contact_cutoff,
            polar_contact_cutoff=(
                polar_contact_cutoff
            ),
        )

    second_level_matrices = (
        assign_second_level_subfamilies(
            all_pose_records=all_pose_records,
            distance_cutoff=(
                subfamily_distance_cutoff
            ),
            orientation_weight=orientation_weight,
            contact_weight=contact_weight,
            minimum_axis_anisotropy=(
                minimum_axis_anisotropy
            ),
            linkage_method=linkage_method,
        )
    )

    logger.info(
        "Assigned orientation/contact subfamilies inside "
        "%d global spatial families",
        len(
            second_level_matrices
        ),
    )

    # --------------------------------------------------
    # Flat pose table
    # --------------------------------------------------

    pose_rows = []

    for record in all_pose_records:
        centroid = record[
            "centroid"
        ]

        pose_rows.append(
            {
                "source_name": record[
                    "source_name"
                ],
                "pose_number": record[
                    "pose_number"
                ],
                "local_cluster_id": record[
                    "local_cluster_id"
                ],
                "global_family_id": record[
                    "global_family_id"
                ],
                "global_spatial_family_id": record[
                    "global_spatial_family_id"
                ],
                "orientation_contact_subfamily_id": record[
                    "orientation_contact_subfamily_id"
                ],
                "is_family_medoid": record[
                    "is_family_medoid"
                ],
                "centroid_x": float(
                    centroid[
                        0
                    ]
                ),
                "centroid_y": float(
                    centroid[
                        1
                    ]
                ),
                "centroid_z": float(
                    centroid[
                        2
                    ]
                ),
                "principal_axis_x": float(
                    record[
                        "principal_axis"
                    ][
                        0
                    ]
                ),
                "principal_axis_y": float(
                    record[
                        "principal_axis"
                    ][
                        1
                    ]
                ),
                "principal_axis_z": float(
                    record[
                        "principal_axis"
                    ][
                        2
                    ]
                ),
                "axis_anisotropy": float(
                    record[
                        "axis_anisotropy"
                    ]
                ),
                "radius_of_gyration": float(
                    record[
                        "radius_of_gyration"
                    ]
                ),
                "n_contact_features": int(
                    len(
                        record[
                            "contact_fingerprint"
                        ]
                    )
                ),
                "contact_fingerprint": ";".join(
                    sorted(
                        record[
                            "contact_fingerprint"
                        ]
                    )
                ),                
                "n_heavy_atoms": int(
                    len(
                        record[
                            "coordinates"
                        ]
                    )
                ),
                "reference_rmsd": record[
                    "reference_rmsd"
                ],
                "density_score": record[
                    "density_score"
                ],
                "docking_score": record[
                    "docking_score"
                ],
            }
        )

    pose_table = pd.DataFrame(
        pose_rows
    )

    pose_table.to_csv(
        output_dir
        / "campaign_aligned_pose_table.csv",
        index=False,
    )

    # --------------------------------------------------
    # Posterior ensemble-support analysis
    # --------------------------------------------------

    ensemble_support_results = None

    if ensemble_support_enabled:
        logger.info(
            "Running posterior ensemble-support analysis "
            "from objectively inferred spatial and "
            "orientation/contact modes"
        )

        ensemble_support_results = (
            run_ensemble_support_analysis(
                pose_table=pose_table,
                output_dir=(
                    output_dir
                    / "ensemble_support"
                ),
                config=ensemble_support_cfg,
            )
        )

        logger.info(
            "Posterior ensemble-support analysis completed: %s",
            output_dir
            / "ensemble_support",
        )

    pocket_evidence_results = None

    if pocket_evidence_enabled:
        if ensemble_support_results is None:
            raise RuntimeError(
                "Pocket-evidence analysis requires the "
                "ensemble-support analysis."
            )

        joint_adjusted_summary = (
            ensemble_support_results.get(
                "_joint_adjusted_summary_df"
            )
        )

        if joint_adjusted_summary is None:
            raise RuntimeError(
                "Ensemble-support analysis did not return the "
                "joint adjusted posterior summary."
            )

        # Export one reference receptor PDB into the pocket-evidence
        # directory so its PML remains portable.
        pocket_output_dir = (
            output_dir
            / "pocket_evidence"
        )

        pocket_output_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        # ----------------------------------------------
        # Retrieve the receptor-only aligned copy of the
        # campaign reference structure.
        # ----------------------------------------------

        reference_receptor_for_pocket = None

        for receptor_record in aligned_receptor_records:
            if (
                receptor_record[
                    "source_name"
                ]
                == reference_name
            ):
                reference_receptor_for_pocket = (
                    receptor_record[
                        "structure"
                    ]
                )

                break

        if reference_receptor_for_pocket is None:
            raise RuntimeError(
                "Could not locate the aligned receptor-only "
                f"reference structure '{reference_name}' for "
                "pocket-evidence analysis."
            )

        receptor_pdb = (
            pocket_output_dir
            / "reference_receptor.pdb"
        )

        with structure.StructureWriter(
            str(
                receptor_pdb
            )
        ) as writer:
            writer.append(
                reference_receptor_for_pocket
            )

        pocket_evidence_results = (
            run_pocket_evidence_analysis(
                pose_table=pose_table,
                joint_summary=(
                    joint_adjusted_summary
                ),
                reference_receptor=(
                    reference_receptor_for_pocket
                ),
                reference_chain_id=(
                    campaign_reference_chain_id
                ),
                reference_receptor_pdb_name=(
                    receptor_pdb.name
                ),
                output_dir=(
                    pocket_output_dir
                ),
                config=pocket_evidence_cfg,
            )
        )

        logger.info(
            "Generated posterior-weighted pocket evidence "
            "analysis: %s",
            pocket_output_dir,
        )

    contact_summary = (
        second_level_contact_summary(
            pose_table
        )
    )

    contact_summary.to_csv(
        output_dir
        / "second_level_contact_summary.csv",
        index=False,
    )    

    second_level_composition = (
        pose_table
        .groupby(
            [
                "global_spatial_family_id",
                "orientation_contact_subfamily_id",
                "source_name",
            ],
            dropna=False,
        )
        .size()
        .rename(
            "n_poses"
        )
        .reset_index()
    )

    save_second_level_family_heatmap(
        pose_table=pose_table,
        output_file=(
            output_dir
            / "second_level_family_heatmap.png"
        ),
    )    

    second_level_totals = (
        second_level_composition
        .groupby(
            [
                "global_spatial_family_id",
                "source_name",
            ]
        )[
            "n_poses"
        ]
        .transform(
            "sum"
        )
    )

    second_level_composition[
        "fraction_within_spatial_family"
    ] = (
        second_level_composition[
            "n_poses"
        ]
        / second_level_totals
    )

    second_level_composition.to_csv(
        output_dir
        / "second_level_subfamily_composition.csv",
        index=False,
    )

    second_level_support = (
        pose_table
        .groupby(
            [
                "global_spatial_family_id",
                "orientation_contact_subfamily_id",
            ],
            dropna=False,
        )
        .agg(
            n_poses=(
                "pose_number",
                "size",
            ),
            n_supporting_structures=(
                "source_name",
                "nunique",
            ),
            mean_axis_anisotropy=(
                "axis_anisotropy",
                "mean",
            ),
            mean_contact_features=(
                "n_contact_features",
                "mean",
            ),
        )
        .reset_index()
        .sort_values(
            [
                "n_supporting_structures",
                "n_poses",
            ],
            ascending=[
                False,
                False,
            ],
        )
    )

    second_level_support.to_csv(
        output_dir
        / "second_level_subfamily_support.csv",
        index=False,
    )

    # --------------------------------------------------
    # Family composition and support
    # --------------------------------------------------

    composition = (
        structure_family_composition(
            pose_table
        )
    )

    composition.to_csv(
        output_dir
        / "structure_global_family_composition.csv",
        index=False,
    )

    ambiguity = (
        structure_ambiguity_summary(
            composition
        )
    )

    ambiguity.to_csv(
        output_dir
        / "structure_global_family_ambiguity.csv",
        index=False,
    )

    family_support = (
        family_support_summary(
            pose_table
        )
    )

    family_support.to_csv(
        output_dir
        / "global_family_support.csv",
        index=False,
    )

    save_structure_family_heatmap(
        composition=composition,
        ambiguity=ambiguity,
        output_file=(
            output_dir
            / "structure_global_family_heatmap.png"
        ),
    )

    # --------------------------------------------------
    # Target, structure and family probability densities
    #
    # Construct one normalized mean density per structure first.
    # The campaign consensus then gives every structure equal total
    # weight, irrespective of whether GlideEM returned 6 or 93 poses.
    # --------------------------------------------------

    structure_densities = {}
    structure_pose_variabilities = {}

    source_names = sorted(
        pose_table[
            "source_name"
        ].unique()
    )

    for source_name in source_names:
        source_pose_densities = [
            record[
                "density"
            ]
            for record in all_pose_records
            if record[
                "source_name"
            ] == source_name
        ]

        if not source_pose_densities:
            raise RuntimeError(
                f"No pose densities were found for "
                f"{source_name}."
            )

        # Mean probability field for this structure.
        structure_densities[
            source_name
        ] = pooled_probability_density(
            source_pose_densities
        )

        # Voxelwise spread among this structure's generated poses.
        structure_pose_variabilities[
            source_name
        ] = density_variability(
            source_pose_densities
        )

    # Equal-structure campaign consensus.
    target_density = pooled_probability_density(
        list(
            structure_densities.values()
        )
    )

    # Voxelwise variation between structure-level mean landscapes.
    between_structure_variability = density_variability(
        list(
            structure_densities.values()
        )
    )

    # Typical voxelwise variation within an individual structure.
    mean_within_structure_variability = np.mean(
        np.stack(
            list(
                structure_pose_variabilities.values()
            ),
            axis=0,
        ),
        axis=0,
    ).astype(
        np.float32
    )

    save_numpy_density(
        target_density,
        output_dir
        / "target_consensus_probability.npy",
    )

    save_numpy_density(
        between_structure_variability,
        output_dir
        / "between_structure_variability.npy",
    )

    save_numpy_density(
        mean_within_structure_variability,
        output_dir
        / "mean_within_structure_variability.npy",
    )

    logger.info(
        "Generated equal-structure campaign density from "
        "%d structure-level probability fields",
        len(
            structure_densities
        ),
    )

    family_densities = {}

    for source_name in source_names:
        source_density = structure_densities[
            source_name
        ]

        save_numpy_density(
            source_density,
            output_dir
            / (
                f"structure_{source_name}"
                "_probability.npy"
            ),
        )

        if export_cfg.get(
            "difference_maps",
            True,
        ):
            difference = density_difference(
                source_density,
                target_density,
            )

            save_numpy_density(
                difference,
                output_dir
                / (
                    f"structure_{source_name}"
                    "_minus_target.npy"
                ),
            )

        # if export_cfg.get(
        #     "difference_maps",
        #     True,
        # ):
        #     difference = density_difference(
        #         source_density,
        #         target_density,
        #     )

        #     save_numpy_density(
        #         difference,
        #         output_dir
        #         / (
        #             f"structure_{source_name}"
        #             "_minus_target.npy"
        #         ),
        #     )

        # if export_cfg.get(
        #     "per_structure_maps",
        #     True,
        # ):
        #     try:
        #         write_ccp4_map(
        #             source_density,
        #             common_grid,
        #             output_dir
        #             / (
        #                 f"structure_{source_name}"
        #                 "_probability.ccp4"
        #             ),
        #         )

        #     except Exception as exc:
        #         logger.warning(
        #             "[%s] Could not export structure CCP4 map: %s",
        #             source_name,
        #             exc,
        #         )

    for family_id in sorted(
        pose_table[
            "global_family_id"
        ].unique()
    ):
        family_indices = [
            index
            for index, record in enumerate(
                all_pose_records
            )
            if record[
                "global_family_id"
            ] == family_id
        ]

        family_density = (
            pooled_probability_density(
                [
                    pose_densities[
                        index
                    ]
                    for index in family_indices
                ]
            )
        )

        family_densities[
            int(
                family_id
            )
        ] = family_density

        save_numpy_density(
            family_density,
            output_dir
            / (
                f"global_family_{int(family_id):03d}"
                "_probability.npy"
            ),
        )

        # if export_cfg.get(
        #     "per_family_maps",
        #     True,
        # ):
        #     try:
        #         write_ccp4_map(
        #             family_density,
        #             common_grid,
        #             output_dir
        #             / (
        #                 f"global_family_{int(family_id):03d}"
        #                 "_probability.ccp4"
        #             ),
        #         )

        #     except Exception as exc:
        #         logger.warning(
        #             "[Family %d] Could not export CCP4 map: %s",
        #             int(
        #                 family_id
        #             ),
        #             exc,
        #         )

    # --------------------------------------------------
    # Between-structure Jensen-Shannon divergence
    # --------------------------------------------------

    js_matrix = structure_js_distance_matrix(
        structure_densities
    )

    js_matrix.to_csv(
        output_dir
        / "structure_landscape_js_distances.csv"
    )

    save_structure_js_heatmap(
        js_matrix=js_matrix,
        output_file=(
            output_dir
            / "structure_landscape_js_heatmap.png"
        ),
    )

    # --------------------------------------------------
    # Within vs between structure spatial dispersion
    # --------------------------------------------------

    (
        structure_dispersion,
        overall_dispersion,
    ) = centroid_dispersion_summary(
        pose_table
    )

    structure_dispersion.to_csv(
        output_dir
        / "structure_centroid_dispersion.csv",
        index=False,
    )

    overall_dispersion.to_csv(
        output_dir
        / "campaign_centroid_variance_decomposition.csv",
        index=False,
    )

    # --------------------------------------------------
    # Export all aligned poses
    # --------------------------------------------------

    if export_cfg.get(
        "aligned_poses",
        True,
    ):
        write_aligned_pose_archive(
            pose_records=all_pose_records,
            output_file=(
                output_dir
                / "all_aligned_glideem_poses.maegz"
            ),
        )

    second_level_pymol_script = None

    if export_second_level_pymol_enabled:
        second_level_pymol_script = (
            export_second_level_pymol(
                pose_records=all_pose_records,
                aligned_receptor_file=(
                    aligned_receptor_file
                ),
                alignment_reference_name=(
                    reference_name
                ),
                output_dir=(
                    output_dir
                    / "second_level_pymol"
                ),
            )
        )

        logger.info(
            "Second-level PyMOL landscape exported: %s",
            second_level_pymol_script,
        )

    logger.info(
        "Campaign GlideEM landscape complete: "
        "%d poses, %d structures, %d global families",
        len(
            pose_table
        ),
        pose_table[
            "source_name"
        ].nunique(),
        pose_table[
            "global_family_id"
        ].nunique(),
    )

    return {
        "pose_table": (
            output_dir
            / "campaign_aligned_pose_table.csv"
        ),
        "aligned_receptors": (
            output_dir
            / "all_aligned_receptors.maegz"
        ),
        "family_composition": (
            output_dir
            / "structure_global_family_composition.csv"
        ),
        "family_support": (
            output_dir
            / "global_family_support.csv"
        ),
        "ambiguity": (
            output_dir
            / "structure_global_family_ambiguity.csv"
        ),
        "alignment": (
            output_dir
            / "receptor_pocket_alignment.csv"
        ),
        "target_density": (
            output_dir
            / "target_consensus_probability.npy"
        ),
        "between_structure_variability": (
            output_dir
            / "between_structure_variability.npy"
        ),
        "mean_within_structure_variability": (
            output_dir
            / "mean_within_structure_variability.npy"
        ),
        "js_distances": (
            output_dir
            / "structure_landscape_js_distances.csv"
        ),
        "variance_decomposition": (
            output_dir
            / "campaign_centroid_variance_decomposition.csv"
        ),
        "second_level_composition": (
            output_dir
            / "second_level_subfamily_composition.csv"
        ),
        "second_level_support": (
            output_dir
            / "second_level_subfamily_support.csv"
        ),
        "alternative_pose_comparison": (
            output_dir
            / "alternative_pose_mode_comparison.csv"
        ),
        "second_level_pymol": (
            second_level_pymol_script
        ),        
        "ensemble_support": (
            ensemble_support_results
        ),
        "pocket_evidence": (
            pocket_evidence_results
        ),
    }
    