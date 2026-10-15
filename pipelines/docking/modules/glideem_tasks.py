from __future__ import annotations

import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, pearsonr, spearmanr
import matplotlib.pyplot as plt
from schrodinger import structure

from pipeline.logger import setup_logger
from pipeline.task_registry import register_task

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

AA = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","MSE","HID","HIE","HIP","ASH","GLH","LYN","CYX"}
IGNORE = {"HOH","WAT","SOL","DOD","NA","K","CL","CA","MG","ZN","MN","FE","SO4","PO4","GOL","EDO","PEG"}
ID_FIELDS = ("s_user_compound_id", "Compound ID", "Compound_ID", "ID", "Name", "s_m_title")


def norm(value):
    return re.sub(r"[^A-Z0-9]+", "_", str(value).upper()).strip("_")

def property_value(st, requested):
    target = norm(requested)
    for key, value in st.property.items():
        nk = norm(key)
        if nk == target or nk.endswith("_" + target) or target in nk:
            return value
    return None


def compound_id(st):
    """
    Return a stable compound identifier from a Maestro structure.

    Identifier matching is deliberately strict. In particular, generic
    substring matching must not be used for 'ID', because properties such
    as s_i_glide_gridfile also contain the letters 'id'.
    """

    # Exact Maestro property names in priority order.
    exact_keys = (
        "s_user_compound_id",
        "s_sd_ID",
        "s_sd_Compound ID",
        "s_sd_Compound_ID",
        "s_sd_Name",
        "s_m_title",
    )

    for key in exact_keys:
        value = st.property.get(key)

        if value is not None:
            value = str(value).strip()

            if value:
                return value

    # Fallback using normalized property names, but only exact matches
    # or suffix matches. Never use loose substring matching here.
    wanted = (
        "USER_COMPOUND_ID",
        "COMPOUND_ID",
        "ID",
        "NAME",
        "M_TITLE",
    )

    for wanted_name in wanted:

        for key, value in st.property.items():

            normalized_key = norm(key)

            if (
                normalized_key == wanted_name
                or normalized_key.endswith(
                    f"_{wanted_name}"
                )
            ):
                value = str(value).strip()

                if value:
                    return value

    title = st.title.strip()

    return title or None


def glide_score(st):
    props = {norm(k): v for k, v in st.property.items()}
    preferred = ("R_I_GLIDE_GSCORE", "R_I_DOCKING_SCORE", "R_GLIDE_GSCORE", "R_I_GLIDE_DOCKING_SCORE")
    for key in preferred:
        try:
            if key in props:
                return float(props[key])
        except (TypeError, ValueError):
            pass
    for key, value in props.items():
        if "GLIDE" in key and ("GSCORE" in key or "DOCKING_SCORE" in key):
            try:
                return float(value)
            except (TypeError, ValueError):
                pass
    return None


def is_protein_atom(atom):
    try:
        return bool(atom.isProtein())
    except Exception:
        return (atom.pdbres or "").strip().upper() in AA


def heavy_atoms(st, indices):
    return sum(st.atom[i].atomic_number > 1 for i in indices)

def annotate_ligprep_compound_ids(
    input_sdf: Path,
    ligprep_file: Path,
) -> None:
    """
    Add a durable s_user_compound_id property to LigPrep structures.

    Assumes one LigPrep output structure per input structure.
    The assumption is explicitly validated before writing.
    """

    source_structures = list(
        structure.StructureReader(
            str(input_sdf)
        )
    )

    prepared_structures = list(
        structure.StructureReader(
            str(ligprep_file)
        )
    )

    if len(source_structures) != len(prepared_structures):
        raise RuntimeError(
            "Cannot safely transfer compound IDs from the "
            "input SDF to LigPrep output: "
            f"{len(source_structures)} input structures but "
            f"{len(prepared_structures)} prepared structures. "
            "A provenance-based mapping is required."
        )

    seen_ids = set()

    for index, (source, prepared) in enumerate(
        zip(
            source_structures,
            prepared_structures,
        ),
        start=1,
    ):
        cid = property_value(
            source,
            "ID",
        )

        if cid is None or not str(cid).strip():
            raise RuntimeError(
                f"Input ligand {index} has no "
                f"'ID' property."
            )

        cid = str(cid).strip()

        if cid in seen_ids:
            raise RuntimeError(
                f"Duplicate compound ID in input SDF: "
                f"{cid}"
            )

        seen_ids.add(cid)

        prepared.property[
            "s_user_compound_id"
        ] = cid

    tmp_file = (
        ligprep_file.parent /
        f"{ligprep_file.stem}_annotated.maegz"
    )

    with structure.StructureWriter(
        str(tmp_file)
    ) as writer:

        for prepared in prepared_structures:
            writer.append(prepared)

    tmp_file.replace(
        ligprep_file
    )

    logger.info(
        f"Annotated {len(prepared_structures)} "
        f"LigPrep structures with "
        f"s_user_compound_id"
    )

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

@register_task("prepare_glideem_receptor", category="GlideEM", description="Prepare a receptor using Protein Preparation Wizard")
def prepare_glideem_receptor(backend, ligand, config, **kwargs):
    job = kwargs["job"]
    work = Path(kwargs["work_dir"])
    source = Path(job["receptor"]).expanduser().resolve()
    output = work / f"{job['job_name']}_prepared.maegz"
    if output.exists() and not config.get("force", False):
        return output
    prep = config["glideem"].get("prepare_receptor", {})
    cmd = [backend.exe("prepwizard"), source, output, "-WAIT", "-HOST", backend.host,
           "-fillsidechains", "-disulfides", "-rehtreat", "-propka_pH", str(prep.get("pH", 7.4)),
           "-rmsd", str(prep.get("rmsd", 0.3))]
    backend.run(cmd, work, work / f"{job['job_name']}_prepwizard.log")
    if not output.exists():
        raise RuntimeError(f"PrepWizard did not create {output}")
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
    ligand_asl = (
        f"res.ptype {ligand_resname}"
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

    logger.info(
        f"Starting GlideEM for {job['job_name']}: "
        f"ligand_asl='{ligand_asl}', "
        f"mode={mode}, nposes={nposes}"
    )

    backend.run_glideem(
        receptor=receptor,
        density_map=density_map,
        ligand_asl=ligand_asl,
        job_name=job_name,
        cwd=work,
        mode=mode,
        nposes=nposes,
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


def extract_complexes(
    archive: Path,
    out_root: Path,
    source_name: str,
    ligand_resname: str,
    force: bool = False
):

    ligand_resname = (
        ligand_resname
        .strip()
        .upper()
    )

    jobs = []
    pose = 0
    for ct_number, st in enumerate(structure.StructureReader(str(archive)), 1):
        protein = [a.index for a in st.atom if is_protein_atom(a)]
        if not protein:
            logger.warning("%s CT %d has no protein; skipped", archive.name, ct_number)
            continue
        residues = defaultdict(list)
        for atom in st.atom:
            if is_protein_atom(atom):
                continue
            rn = (atom.pdbres or "").strip().upper()
            if rn not in IGNORE:
                residues[((atom.chain or "").strip(), int(atom.resnum), rn)].append(atom.index)
        
        candidates = [
            (heavy_atoms(st, idx), key, idx)
            for key, idx in residues.items()
            if (
                key[2] == ligand_resname
                and heavy_atoms(st, idx) >= 6
            )
        ]

        if not candidates:
            logger.warning(
                "%s CT %d has no '%s' ligand; skipped",
                archive.name,
                ct_number,
                ligand_resname,
            )
            continue

        pose += 1
        name = f"{source_name}_pose{pose:02d}"
        directory = out_root / name
        directory.mkdir(parents=True, exist_ok=True)
        receptor = directory / f"{name}_receptor.maegz"
        reference = directory / f"{name}_reference_ligand.maegz"
        complex_file = directory / f"{name}_complex.maegz"
        grid = directory / f"{name}_grid.zip"
        metadata = directory / "complex_metadata.json"

        nheavy, residue, ligand_indices = max(candidates, key=lambda item: item[0])

        logger.info(
            "%s CT %d: found '%s' ligand "
            "(chain=%s, resnum=%d, heavy_atoms=%d) -> pose %02d",
            archive.name,
            ct_number,
            ligand_resname,
            residue[0] or "-",
            residue[1],
            nheavy,
            pose,
        )

        xyz = np.array([[st.atom[i].x, st.atom[i].y, st.atom[i].z] for i in ligand_indices], dtype=float)
        if force or not all(p.exists() for p in (receptor, reference, complex_file, metadata)):
            lig = st.extract(ligand_indices)
            lig.title = name
            lig.property["s_user_reference_complex"] = name
            with structure.StructureWriter(str(receptor)) as writer:
                writer.append(st.extract(protein))
            with structure.StructureWriter(str(reference)) as writer:
                writer.append(lig)
            with structure.StructureWriter(str(complex_file)) as writer:
                writer.append(st)
            metadata.write_text(json.dumps({
                "job_name": name, "source_archive": str(archive.resolve()), "source_CT": ct_number,
                "ligand_chain": residue[0], "ligand_resnum": residue[1], "ligand_resname": residue[2],
                "ligand_heavy_atoms": nheavy, "grid_center": xyz.mean(axis=0).tolist()
            }, indent=2))
        meta = json.loads(metadata.read_text())
        jobs.append({"name": name, "dir": directory, "receptor": receptor, "reference": reference,
                     "complex": complex_file, "grid": grid, "center": meta["grid_center"]})
    return jobs


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


def analyse_correlations(complex_jobs, poseviewers, ligand_file, activity_column, output_dir, min_pairs=5):
    activities = {}
    for st in structure.StructureReader(str(ligand_file)):
        try:
            value = float(property_value(st, activity_column))
        except (TypeError, ValueError):
            continue
        if math.isfinite(value):
            activities.setdefault(compound_id(st), value)
    if not activities:
        raise RuntimeError(f"No numeric values found for activity property '{activity_column}'")

    by_name = {j["name"]: j for j in complex_jobs}
    details = []
    for name, pv in poseviewers.items():
        best = {}
        for st in structure.StructureReader(str(pv)):
            cid = compound_id(st)
            score = glide_score(st)
            if cid in activities and score is not None and (cid not in best or score < best[cid]):
                best[cid] = score
        for cid, score in best.items():
            details.append({"complex_pose": name, "compound_id": cid, "glide_score": score,
                            "activity": activities[cid], "activity_column": activity_column,
                            "poseviewer_path": str(Path(pv).resolve()),
                            "reference_ligand_path": str(by_name[name]["reference"].resolve()),
                            "grid_path": str(by_name[name]["grid"].resolve())})
    activity_ids = set(
        activities
    )

    docked_ids = set()
    scored_ids = set()
    matched_ids = set()

    by_name = {
        j["name"]: j
        for j in complex_jobs
    }

    details = []

    docked_ids = set()
    scored_ids = set()
    matched_ids = set()

    for name, pv in poseviewers.items():

        best = {}

        for st in structure.StructureReader(
            str(pv)
        ):
            cid = compound_id(st)
            score = glide_score(st)

            if cid:
                docked_ids.add(cid)

            if cid and score is not None:
                scored_ids.add(cid)

            if (
                cid in activities
                and score is not None
                and (
                    cid not in best
                    or score < best[cid]
                )
            ):
                best[cid] = score
                matched_ids.add(cid)

        logger.info(
            "[%s] Docked IDs with scores: %d; "
            "matched to activity: %d",
            name,
            len(best),
            sum(
                cid in activities
                for cid in best
            ),
        )

        for cid, score in best.items():

            details.append(
                {
                    "complex_pose": name,
                    "compound_id": cid,
                    "glide_score": score,
                    "activity": activities[cid],
                    "activity_column": activity_column,
                    "poseviewer_path": str(
                        Path(pv).resolve()
                    ),
                    "reference_ligand_path": str(
                        by_name[name][
                            "reference"
                        ].resolve()
                    ),
                    "grid_path": str(
                        by_name[name][
                            "grid"
                        ].resolve()
                    ),
                }
            )

    if not details:
        raise RuntimeError(
            "No matched docking score/activity pairs. "
            f"Activity IDs={len(activities)}, "
            f"docked IDs={len(docked_ids)}, "
            f"scored IDs={len(scored_ids)}, "
            f"matched IDs={len(matched_ids)}."
        )
    
    if not details:
        raise RuntimeError("No matched docking score/activity pairs")

    output_dir.mkdir(parents=True, exist_ok=True)

    detail_df = pd.DataFrame(details)
    detail_df.to_csv(output_dir / "matched_scores_and_activities.csv", index=False)
    summaries = []
    for name, group in detail_df.groupby("complex_pose"):
        if len(group) < min_pairs:
            logger.warning("%s has only %d matched pairs; excluded", name, len(group))
            continue
        x = group["glide_score"].to_numpy(float)
        y = group["activity"].to_numpy(float)
        p, s, k = pearsonr(x, y), spearmanr(x, y), kendalltau(x, y)
        # pActivity should normally anticorrelate with docking score, hence larger -(rP+rS) is better.
        summaries.append({"complex_pose": name, "n": len(group), "pearson_r": p.statistic, "pearson_p": p.pvalue,
                          "spearman_rho": s.statistic, "spearman_p": s.pvalue, "kendall_tau": k.statistic,
                          "kendall_p": k.pvalue, "rank_score": -(p.statistic + s.statistic)})
    summary = pd.DataFrame(summaries)
    if summary.empty:
        raise RuntimeError("No complex had enough matched pairs for correlation analysis")
    summary = summary.sort_values(["rank_score", "spearman_rho", "pearson_r"], ascending=[False, True, True]).reset_index(drop=True)
    summary.insert(0, "rank", np.arange(1, len(summary) + 1))
    summary.to_csv(output_dir / "complex_pose_correlations.csv", index=False)

    for name, group in detail_df.groupby("complex_pose"):

        if len(group) < min_pairs:
            logger.warning(
                "%s has only %d matched pairs; "
                "correlation plot skipped",
                name,
                len(group),
            )
            continue

        # Conventional experimental vs predicted orientation:
        # experimental values on X, predicted values on Y.
        x = group[
            "activity"
        ].to_numpy(float)

        y = group[
            "glide_score"
        ].to_numpy(float)

        # ---------------------------------------------
        # Correlation statistics
        # ---------------------------------------------

        pearson = pearsonr(x, y)
        spearman = spearmanr(x, y)
        kendall = kendalltau(x, y)

        # ---------------------------------------------
        # Linear regression / best-fit line
        # ---------------------------------------------

        slope, intercept = np.polyfit(
            x,
            y,
            1,
        )

        x_fit = np.linspace(
            x.min(),
            x.max(),
            200,
        )

        y_fit = (
            slope * x_fit
            + intercept
        )

        # ---------------------------------------------
        # Plot
        # ---------------------------------------------

        fig, ax = plt.subplots(
            figsize=(7, 6)
        )

        ax.scatter(
            x,
            y,
            alpha=0.8,
            label="Compounds",
        )

        ax.plot(
            x_fit,
            y_fit,
            linewidth=2,
            label="Linear fit",
        )

        ax.set_title(name)

        ax.set_xlabel(
            activity_column
        )

        ax.set_ylabel(
            "Glide score"
        )

        # ---------------------------------------------
        # Statistics annotation
        # ---------------------------------------------

        stats_text = (
            f"N = {len(group)}\n"
            f"Pearson r = {pearson.statistic:.3f} "
            f"(p = {pearson.pvalue:.2g})\n"
            f"Spearman ρ = {spearman.statistic:.3f} "
            f"(p = {spearman.pvalue:.2g})\n"
            f"Kendall τ = {kendall.statistic:.3f} "
            f"(p = {kendall.pvalue:.2g})"
        )

        ax.text(
            0.03,
            0.97,
            stats_text,
            transform=ax.transAxes,
            verticalalignment="top",
            horizontalalignment="left",
            bbox={
                "boxstyle": "round",
                "facecolor": "white",
                "alpha": 0.85,
            },
        )

        ax.legend()

        fig.tight_layout()

        fig.savefig(
            output_dir /
            f"{name}_score_vs_activity.png",
            dpi=200,
        )

        plt.close(fig)
    best = summary.iloc[0]
    (output_dir / "best_binding_pose.txt").write_text(
        f"Best complex pose: {best.complex_pose}\nN: {best.n}\nPearson r: {best.pearson_r:.6f}\n"
        f"Spearman rho: {best.spearman_rho:.6f}\nKendall tau: {best.kendall_tau:.6f}\n"
    )
    return summary
