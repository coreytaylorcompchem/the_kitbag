import re
import json
import shutil

import numpy as np
import pandas as pd
from collections import defaultdict
from pathlib import Path

from schrodinger import structure

from modules.utils.constants import AA, IGNORE, element_colour_commands

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def is_protein_atom(atom):
    try:
        return bool(atom.isProtein())
    except Exception:
        return (atom.pdbres or "").strip().upper() in AA

def norm(value):
    return re.sub(r"[^A-Z0-9]+", "_", str(value).upper()).strip("_")

def property_value(st, requested):
    target = norm(requested)
    for key, value in st.property.items():
        nk = norm(key)
        if nk == target or nk.endswith("_" + target) or target in nk:
            return value
    return None

def heavy_atoms(st, indices):
    return sum(st.atom[i].atomic_number > 1 for i in indices)

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

def smiles_value(st):
    """
    Retrieve SMILES from a docked Maestro structure.
    """

    preferred = (
        "s_sd_SMILES",
        "s_m_smiles",
        "SMILES",
    )

    for key in preferred:
        value = st.property.get(key)

        if value is not None and str(value).strip():
            return str(value).strip()

    # Fallback for variant property naming
    for key, value in st.property.items():
        if "SMILES" in key.upper():
            if value is not None and str(value).strip():
                return str(value).strip()

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

    # Fallback using normalised property names, but only exact matches
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

            normalised_key = norm(key)

            if (
                normalised_key == wanted_name
                or normalised_key.endswith(
                    f"_{wanted_name}"
                )
            ):
                value = str(value).strip()

                if value:
                    return value

    title = st.title.strip()

    return title or None

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

def annotate_ligprep_csv_metadata(
    ligprep_file: Path,
    metadata_by_id: dict,
    activity_column: str,
) -> None:
    """
    Restore compound ID, SMILES and experimental activity
    onto LigPrep output structures produced from CSV/SMILES input.
    """

    prepared_structures = list(
        structure.StructureReader(
            str(ligprep_file)
        )
    )

    if not prepared_structures:
        raise RuntimeError(
            "LigPrep produced no structures: "
            f"{ligprep_file}"
        )

    annotated = 0
    unmatched = []

    for prepared in prepared_structures:

        cid = (
            prepared.title
            or ""
        ).strip()

        if not cid:
            unmatched.append(
                "<missing title>"
            )
            continue

        metadata = metadata_by_id.get(
            cid
        )

        if metadata is None:
            unmatched.append(
                cid
            )
            continue

        prepared.property[
            "s_user_compound_id"
        ] = cid

        prepared.property[
            "s_user_input_smiles"
        ] = metadata[
            "smiles"
        ]

        activity = metadata[
            "activity"
        ]

        if activity is not None:
            prepared.property[
                f"r_user_{norm(activity_column).lower()}"
            ] = float(
                activity
            )

        annotated += 1

    if unmatched:
        examples = ", ".join(
            unmatched[:10]
        )

        raise RuntimeError(
            "Could not map some LigPrep output structures "
            "back to the activity CSV using their titles. "
            f"Examples: {examples}"
        )

    tmp_file = (
        ligprep_file.parent
        / (
            f"{ligprep_file.stem}"
            "_annotated.maegz"
        )
    )

    with structure.StructureWriter(
        str(tmp_file)
    ) as writer:
        for prepared in (
            prepared_structures
        ):
            writer.append(
                prepared
            )

    tmp_file.replace(
        ligprep_file
    )

    logger.info(
        "Annotated %d LigPrep structures with "
        "CSV compound metadata",
        annotated,
    )

def get_reference_ligand_atoms(
    st,
    ligand_resname: str,
):
    """
    Return atom indices for the reference ligand residue.

    Fails if no matching residue is found or if the residue name
    identifies multiple separate residues.
    """

    ligand_resname = (
        ligand_resname
        .strip()
        .upper()
    )

    residues = defaultdict(list)

    for atom in st.atom:

        resname = (
            atom.pdbres or ""
        ).strip().upper()

        if resname != ligand_resname:
            continue

        key = (
            (atom.chain or "").strip(),
            int(atom.resnum),
            resname,
        )

        residues[key].append(
            atom.index
        )

    if not residues:
        raise RuntimeError(
            f"Reference ligand '{ligand_resname}' "
            "was not found in the input structure."
        )

    if len(residues) > 1:
        locations = ", ".join(
            f"chain={chain or '-'} resnum={resnum}"
            for chain, resnum, _ in residues
        )

        raise RuntimeError(
            f"Reference ligand '{ligand_resname}' occurs "
            f"multiple times: {locations}. "
            "The manifest needs chain/residue-number selectors."
        )

    residue_key, atom_indices = next(
        iter(residues.items())
    )

    return residue_key, atom_indices

def get_reference_ligand_center(
    structure_path: Path,
    ligand_resname: str,
) -> tuple[float, float, float]:
    """
    Calculate the geometric center of the reference ligand
    from its experimental coordinates.
    """

    entries = list(
        structure.StructureReader(
            str(structure_path)
        )
    )

    if len(entries) != 1:
        raise RuntimeError(
            f"Expected one structure in {structure_path}, "
            f"found {len(entries)}"
        )

    st = entries[0]

    _, ligand_atoms = get_reference_ligand_atoms(
        st,
        ligand_resname,
    )

    heavy_indices = [
        index
        for index in ligand_atoms
        if st.atom[index].atomic_number > 1
    ]

    if not heavy_indices:
        raise RuntimeError(
            f"Reference ligand '{ligand_resname}' "
            "contains no heavy atoms."
        )

    xyz = np.array(
        [
            [
                st.atom[index].x,
                st.atom[index].y,
                st.atom[index].z,
            ]
            for index in heavy_indices
        ],
        dtype=float,
    )

    center = xyz.mean(
        axis=0
    )

    return (
        float(center[0]),
        float(center[1]),
        float(center[2]),
    )

def find_glideem_ligand_indices(
    st,
    ligand_resname: str,
    ligand_from_smiles: bool,
):
    """
    Locate the docked ligand in one GlideEM output structure.

    For the original-complex route, require the specified residue name.
    For the SMILES route, use the largest substantial non-protein,
    non-solvent residue.
    """

    ligand_resname = (
        ligand_resname
        .strip()
        .upper()
    )

    residues = defaultdict(list)

    for atom in st.atom:

        if is_protein_atom(atom):
            continue

        resname = (
            atom.pdbres or ""
        ).strip().upper()

        if resname in IGNORE:
            continue

        key = (
            (atom.chain or "").strip(),
            int(atom.resnum),
            resname,
        )

        residues[key].append(
            atom.index
        )

    if ligand_from_smiles:
        candidates = [
            (
                heavy_atoms(st, indices),
                key,
                indices,
            )
            for key, indices in residues.items()
            if heavy_atoms(st, indices) >= 6
        ]

    else:
        candidates = [
            (
                heavy_atoms(st, indices),
                key,
                indices,
            )
            for key, indices in residues.items()
            if (
                key[2] == ligand_resname
                and heavy_atoms(st, indices) >= 6
            )
        ]

    if not candidates:
        return None

    _, _, indices = max(
        candidates,
        key=lambda item: item[0],
    )

    return indices

def ligand_heavy_atom_coordinates(
    st,
    atom_indices,
):
    """
    Return heavy-atom coordinates and element sequence in structure order.
    """

    heavy_indices = [
        index
        for index in atom_indices
        if st.atom[index].atomic_number > 1
    ]

    coordinates = np.array(
        [
            [
                st.atom[index].x,
                st.atom[index].y,
                st.atom[index].z,
            ]
            for index in heavy_indices
        ],
        dtype=float,
    )

    elements = tuple(
        st.atom[index].atomic_number
        for index in heavy_indices
    )

    return coordinates, elements

def coordinate_frame_rmsd(
    coordinates_a,
    coordinates_b,
):
    """
    RMSD without superposition.

    This retains translational, rotational and conformational
    differences within the common receptor coordinate frame.
    """

    if coordinates_a.shape != coordinates_b.shape:
        raise ValueError(
            "Cannot compare poses with different coordinate shapes: "
            f"{coordinates_a.shape} versus {coordinates_b.shape}"
        )

    displacement = (
        coordinates_a
        - coordinates_b
    )

    return float(
        np.sqrt(
            np.mean(
                np.sum(
                    displacement * displacement,
                    axis=1,
                )
            )
        )
    )

def glide_score_components(st):
    """
    Extract useful Glide score components from a docked structure.
    """

    fields = {
        "glide_gscore": "r_i_glide_gscore",
        "docking_score": "r_i_docking_score",
        "glide_emodel": "r_i_glide_emodel",
        "glide_energy": "r_i_glide_energy",
        "glide_evdw": "r_i_glide_evdw",
        "glide_ecoul": "r_i_glide_ecoul",
        "glide_lipo": "r_i_glide_lipo",
        "glide_hbond": "r_i_glide_hbond",
        "glide_erotb": "r_i_glide_erotb",
        "glide_esite": "r_i_glide_esite",
        "glide_einternal": "r_i_glide_einternal",
    }

    result = {}

    for output_name, property_name in fields.items():

        value = st.property.get(
            property_name
        )

        try:
            result[output_name] = (
                float(value)
                if value is not None
                else None
            )

        except (TypeError, ValueError):
            result[output_name] = None

    return result

def glide_smiles_by_lignum(poseviewer: Path):
    """
    Read the Glide CSV corresponding to a PoseViewer and map
    Glide ligand number -> SMILES.
    """

    poseviewer = Path(poseviewer)

    # Example:
    # structure_51052_pose01_mcss_pv.maegz
    # ->
    # structure_51052_pose01_mcss.csv

    suffix = "_pv.maegz"

    if not poseviewer.name.endswith(suffix):
        return {}

    csv_file = poseviewer.with_name(
        poseviewer.name[:-len(suffix)] + ".csv"
    )

    if not csv_file.is_file():
        logger.warning(
            "Glide CSV not found for PoseViewer: %s",
            csv_file,
        )
        return {}

    df = pd.read_csv(
        csv_file,
        usecols=[
            "SMILES",
            "i_i_glide_lignum",
        ],
    )

    result = {}

    for _, row in df.iterrows():

        lignum = row["i_i_glide_lignum"]
        smiles = row["SMILES"]

        if pd.isna(lignum) or pd.isna(smiles):
            continue

        result[int(lignum)] = str(smiles)

    return result

def extract_complexes(
    archive: Path,
    out_root: Path,
    source_name: str,
    ligand_resname: str,
    force: bool = False,
    ligand_from_smiles: bool = False,
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

        # First collect ALL non-protein residues in this CT.
        for atom in st.atom:

            if is_protein_atom(atom):
                continue

            rn = (
                atom.pdbres or ""
            ).strip().upper()

            if rn not in IGNORE:
                residues[
                    (
                        (atom.chain or "").strip(),
                        int(atom.resnum),
                        rn,
                    )
                ].append(
                    atom.index
                )

        # Only after the complete residue dictionary has been
        # constructed should we look for the GlideEM ligand.
        if ligand_from_smiles:

            # GlideEM prepared the ligand from SMILES, so its
            # output residue identity does not need to match
            # the reference ligand in the original PDB.
            candidates = [
                (
                    heavy_atoms(st, indices),
                    key,
                    indices,
                )
                for key, indices in residues.items()
                if heavy_atoms(st, indices) >= 6
            ]

        else:

            # Traditional route: GlideEM used the ligand already
            # present in the prepared receptor.
            candidates = [
                (
                    heavy_atoms(st, indices),
                    key,
                    indices,
                )
                for key, indices in residues.items()
                if (
                    key[2] == ligand_resname
                    and heavy_atoms(st, indices) >= 6
                )
            ]

        if not candidates:

            if ligand_from_smiles:
                logger.warning(
                    "%s CT %d has no substantial "
                    "non-protein ligand; skipped",
                    archive.name,
                    ct_number,
                )

            else:
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

        logger.debug(
            "%s CT %d: found '%s' ligand "
            "(chain=%s, resnum=%d, heavy_atoms=%d) -> pose %02d",
            archive.name,
            ct_number,
            residue[2] or "-",
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
        jobs.append({
            "name": name,
            "source_name": source_name,
            "pose_number": pose,
            "dir": directory,
            "receptor": receptor,
            "reference": reference,
            "complex": complex_file,
            "grid": grid,
            "center": meta["grid_center"],
        })
    return jobs

def export_glideem_pymol_ensemble(
    pose_structures,
    cluster_by_pose,
    cryoem_reference: Path,
    receptor: Path,
    density_map: Path,
    output_dir: Path,
    job_name: str,
):

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    receptor = (
        Path(receptor)
        .expanduser()
        .resolve()
    )

    density_map = (
        Path(density_map)
        .expanduser()
        .resolve()
    )

    if not receptor.is_file():
        raise FileNotFoundError(
            f"Prepared GlideEM receptor not found: "
            f"{receptor}"
        )

    if not density_map.is_file():
        raise FileNotFoundError(
            f"GlideEM density map not found: "
            f"{density_map}"
        )

    receptor_structures = list(
        structure.StructureReader(
            str(receptor)
        )
    )

    if not receptor_structures:
        raise RuntimeError(
            f"No structure found in prepared receptor: "
            f"{receptor}"
        )

    receptor_st = receptor_structures[0]

    receptor_file = (
        output_dir /
        "glideem_receptor.pdb"
    )

    with structure.StructureWriter(
        str(receptor_file)
    ) as writer:
        writer.append(
            receptor_st
        )

    cryoem_st = next(
        structure.StructureReader(
            str(cryoem_reference)
        )
    )

    cryoem_file = (
        output_dir /
        "cryoem_structure.pdb"
    )

    with structure.StructureWriter(
        str(cryoem_file)
    ) as writer:
        writer.append(cryoem_st)

    # Only bundle CCP4 density with the PyMOL session.
    # Full MRC maps can be too large for practical PyMOL use.
    density_file = None

    if density_map.suffix.lower() == ".ccp4":

        density_file = (
            output_dir /
            "glideem_density.ccp4"
        )

        shutil.copy2(
            density_map,
            density_file,
        )

    exported = []

    for record in pose_structures:

        pose_number = record[
            "pose_number"
        ]

        cluster_id = cluster_by_pose[
            pose_number
        ]

        object_name = (
            f"pose{pose_number:02d}"
            f"_cluster{cluster_id:02d}"
        )

        pose_file = (
            output_dir /
            f"{object_name}.pdb"
        )

        with structure.StructureWriter(
            str(pose_file)
        ) as writer:
            writer.append(
                record["structure"]
            )

        exported.append(
            (
                pose_number,
                cluster_id,
                object_name,
                pose_file,
            )
        )

    carbon_colours = [
        "cyan",
        "green",
        "yellow",
        "orange",
        "magenta",
        "marine",
        "salmon",
        "violet",
    ]

    script = [
        "reinitialize",
        "bg_color white",
        "",

        # Prepared protein used for GlideEM.
        "load glideem_receptor.pdb, receptor",
        "hide everything, receptor",
        "show cartoon, receptor",
        "color gray80, receptor",
        "",

        # Original experimental ligand position.
        "load cryoem_structure.pdb, cryoem_structure",
        "show sticks, cryoem_structure",
        *element_colour_commands(
            "cryoem_structure",
            "gray50",
        ),
        "",
    ]

    if density_file is not None:

        script.extend(
            [
                "# Experimental density used for GlideEM.",
                (
                    f"load {density_file.name}, "
                    "cryoem_density"
                ),
                "",
                "# Cryo-EM density surface",
                (
                    "isosurface cryoem_density_surface, "
                    "cryoem_density, 1.0"
                ),
                "color wheat, cryoem_density_surface",
                (
                    "set transparency, 0.70, "
                    "cryoem_density_surface"
                ),
                "",
            ]
        )

    else:
        logger.info(
            "[%s] Density map is not .ccp4 (%s); not bundling it "
            "with the PyMOL session",
            job_name,
            density_map.suffix,
        )

    cluster_objects = defaultdict(list)

    for (
        pose_number,
        cluster_id,
        object_name,
        pose_file,
    ) in exported:

        carbon_colour = carbon_colours[
            (cluster_id - 1)
            % len(carbon_colours)
        ]

        cluster_name = (
            f"cluster{cluster_id:02d}"
        )

        cluster_objects[
            cluster_name
        ].append(
            object_name
        )

        script.extend(
            [
                f"load {pose_file.name}, {object_name}",
                f"show sticks, {object_name}",
                *element_colour_commands(
                    object_name,
                    carbon_colour,
                ),
                "",
            ]
        )

    for cluster_name, objects in sorted(
        cluster_objects.items()
    ):
        script.append(
            f"group {cluster_name}, "
            + " ".join(objects)
        )

    all_pose_objects = [
        object_name
        for objects in cluster_objects.values()
        for object_name in objects
    ]

    if not all_pose_objects:
        raise RuntimeError(
            f"No GlideEM pose objects were exported "
            f"for {job_name}"
        )

    script.extend(
        [
            "",
            "# Selection containing every GlideEM pose.",
            (
                "select glideem_pose_ensemble, "
                + " or ".join(
                    all_pose_objects
                )
            ),
            "",
        ]
    )

    script.extend(
        [
            "# Protein residues surrounding the ligand ensemble.",
            (
                "select binding_site, "
                "byres (receptor within 5.0 of "
                "(cryoem_structure or "
                "glideem_pose_ensemble))"
            ),
            "show sticks, binding_site",
            *element_colour_commands(
                "binding_site",
                "gray80",
            ),
            "",
        ]
    )


    script.extend(
        [
            "# Consistent thin-stick representation.",
            "set stick_radius, 0.15, glideem_pose_ensemble",
            "set stick_radius, 0.15, cryoem_structure",
            "set stick_radius, 0.15, binding_site",
            "",

            "# Initial view centred on the ligand ensemble.",
            "orient glideem_pose_ensemble",
            "zoom glideem_pose_ensemble, 8",
            "",
        ]
    )
    
    script_file = (
        output_dir /
        "load_glideem_ensemble.pml"
    )

    script_file.write_text(
        "\n".join(script)
        + "\n"
    )

    return script_file

def prepare_activity_csv_for_ligprep(
    input_csv: Path,
    work_dir: Path,
    smiles_column: str,
    id_column: str,
    activity_column: str,
) -> tuple[Path, dict]:
    """
    Validate an activity CSV and create a SMILES file for LigPrep.

    Returns:
        smiles_file
        metadata_by_id

    metadata_by_id preserves the original compound ID and
    experimental activity so they can be restored after LigPrep.
    """

    df = pd.read_csv(
        input_csv
    )

    required_columns = {
        smiles_column,
        id_column,
        activity_column,
    }

    missing = (
        required_columns
        - set(df.columns)
    )

    if missing:
        raise ValueError(
            "Activity ligand CSV is missing required "
            f"columns: {sorted(missing)}"
        )

    if df.empty:
        raise ValueError(
            f"Activity ligand CSV is empty: {input_csv}"
        )

    metadata_by_id = {}

    smiles_lines = []

    for row_number, row in df.iterrows():

        cid = row[
            id_column
        ]

        smiles = row[
            smiles_column
        ]

        activity = row[
            activity_column
        ]

        if pd.isna(cid) or not str(
            cid
        ).strip():
            raise ValueError(
                "Missing compound ID in activity CSV "
                f"at row {row_number + 2}"
            )

        if pd.isna(smiles) or not str(
            smiles
        ).strip():
            raise ValueError(
                "Missing SMILES in activity CSV "
                f"at row {row_number + 2}"
            )

        cid = str(
            cid
        ).strip()

        smiles = str(
            smiles
        ).strip()

        if cid in metadata_by_id:
            raise ValueError(
                "Duplicate compound ID in activity CSV: "
                f"{cid}"
            )

        if pd.isna(activity):
            activity_value = None

        else:
            try:
                activity_value = float(
                    activity
                )

            except (TypeError, ValueError):
                raise ValueError(
                    "Non-numeric activity for compound "
                    f"{cid}: {activity!r}"
                )

        metadata_by_id[
            cid
        ] = {
            "activity": activity_value,
            "smiles": smiles,
        }

        # Standard SMILES file:
        #
        # SMILES <whitespace> title
        #
        # The title becomes our stable compound identifier.
        smiles_lines.append(
            f"{smiles}\t{cid}"
        )

    smiles_file = (
        work_dir
        / "activity_ligands_input.smi"
    )

    smiles_file.write_text(
        "\n".join(
            smiles_lines
        )
        + "\n"
    )

    return (
        smiles_file,
        metadata_by_id,
    )