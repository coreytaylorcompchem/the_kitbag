from __future__ import annotations

import csv
import re

import numpy as np
import pandas as pd

import json
from pathlib import Path

from modules.utils.constants import element_colour_commands

from modules.utils.glideem_helpers import (
    compound_id,
)

from modules.utils.glideem_campaign_analysis import (
    transform_structure_copy,
)

from schrodinger import structure


ELEMENT_COLOURS = {
    "N": "blue",
    "O": "red",
    "S": "yellow",
    "P": "orange",
    "H": "white",
    "F": "palecyan",
    "Cl": "green",
    "Br": "brown",
    "I": "purple",
}


def safe_pymol_name(
    value,
    fallback: str = "unknown",
) -> str:
    safe_name = re.sub(
        r"[^A-Za-z0-9_]+",
        "_",
        str(value),
    ).strip("_")

    return safe_name or fallback


def write_structure(
    st,
    output_file: Path,
) -> Path:
    output_file = (
        Path(output_file)
        .expanduser()
        .resolve()
    )

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with structure.StructureWriter(
        str(output_file)
    ) as writer:
        writer.append(
            st
        )

    return output_file


def read_first_structure(
    input_file: Path,
):
    input_file = (
        Path(input_file)
        .expanduser()
        .resolve()
    )

    structures = list(
        structure.StructureReader(
            str(input_file)
        )
    )

    if not structures:
        raise RuntimeError(
            f"No structure found in {input_file}"
        )

    return structures[0]

def write_pymol_script(
    commands: list[str],
    output_file: Path,
) -> Path:
    output_file = (
        Path(output_file)
        .expanduser()
        .resolve()
    )

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_file.write_text(
        "\n".join(commands)
        + "\n"
    )

    return output_file


def export_mcss_pose_bundle(
    receptor,
    ranked_poses,
    cryoem_reference: Path,
    glideem_reference: Path,
    output_dir: Path,
    binding_site_cutoff: float = 5.0,
) -> Path:
    """
    Export an MCSS visualization bundle.

    ranked_poses must contain dictionaries with:
      rank, score, compound_id, structure
    """

    output_dir = (
        Path(output_dir)
        .expanduser()
        .resolve()
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    write_structure(
        receptor,
        output_dir /
        "receptor.pdb",
    )

    cryoem_structure = read_first_structure(
        cryoem_reference
    )

    glideem_structure = read_first_structure(
        glideem_reference
    )

    write_structure(
        cryoem_structure,
        output_dir /
        "cryoem_structure.pdb",
    )

    write_structure(
        glideem_structure,
        output_dir /
        "glideem_reference_pose.pdb",
    )

    exported = []

    for pose in ranked_poses:
        rank = int(
            pose["rank"]
        )

        compound_id = pose[
            "compound_id"
        ]

        safe_id = safe_pymol_name(
            compound_id
        )

        object_name = (
            f"rank{rank:02d}_{safe_id}"
        )

        pose_file = write_structure(
            pose["structure"],
            output_dir /
            f"{object_name}.pdb",
        )

        exported.append(
            {
                **pose,
                "object_name": object_name,
                "pose_file": pose_file,
            }
        )

    ranking_file = (
        output_dir /
        f"top{len(exported)}.csv"
    )

    with ranking_file.open(
        "w",
        newline="",
    ) as handle:
        writer = csv.writer(
            handle
        )

        writer.writerow(
            [
                "rank",
                "compound_id",
                "glide_score",
                "pymol_object",
                "pose_file",
            ]
        )

        for pose in exported:
            writer.writerow(
                [
                    pose["rank"],
                    pose["compound_id"],
                    pose["score"],
                    pose["object_name"],
                    pose["pose_file"].name,
                ]
            )

    commands = [
        "reinitialize",
        "bg_color white",
        "",
        "load receptor.pdb, receptor",
        "hide everything",
        "show cartoon, receptor",
        "color gray80, receptor",
        "",
        "load cryoem_structure.pdb, cryoem_structure",
        "show sticks, cryoem_structure",
        *element_colour_commands(
            "cryoem_structure",
            "marine",
        ),
        "",
        (
            "load glideem_reference_pose.pdb, "
            "glideem_reference_pose"
        ),
        "show sticks, glideem_reference_pose",
        *element_colour_commands(
            "glideem_reference_pose",
            "salmon",
        ),
        "",
    ]

    carbon_colours = [
        "cyan",
        "green",
        "yellow",
        "orange",
        "magenta",
    ]

    ligand_objects = []

    for pose in exported:
        object_name = pose[
            "object_name"
        ]

        pose_file = pose[
            "pose_file"
        ]

        rank = int(
            pose["rank"]
        )

        carbon_colour = carbon_colours[
            (rank - 1)
            % len(carbon_colours)
        ]

        ligand_objects.append(
            object_name
        )

        commands.extend(
            [
                (
                    f"load {pose_file.name}, "
                    f"{object_name}"
                ),
                f"show sticks, {object_name}",
                *element_colour_commands(
                    object_name,
                    carbon_colour,
                ),
                "",
            ]
        )

    if not ligand_objects:
        raise RuntimeError(
            "No MCSS ligand objects were exported."
        )

    # Keep the individual docked ligands inside an expandable
    # top-level group.
    commands.append(
        "group top_mcss_poses, "
        + " ".join(
            ligand_objects
        )
    )

    # Use a separate aggregate selection for display commands.
    # This avoids relying on the group itself as a selection.
    commands.extend(
        [
            "",
            (
                "select mcss_pose_ensemble, "
                + " or ".join(
                    ligand_objects
                )
            ),
            "",
        ]
    )

    ligand_selection = (
        "cryoem_structure or "
        "glideem_reference_pose or "
        "mcss_pose_ensemble"
    )

    commands.extend(
        [
            "",
            "# Protein residues surrounding all displayed ligands.",
            (
                "select binding_site, "
                f"byres (receptor within "
                f"{binding_site_cutoff:.1f} of "
                f"({ligand_selection}))"
            ),
            "show sticks, binding_site",
            *element_colour_commands(
                "binding_site",
                "gray80",
            ),
            "",
            "# Consistent thin-stick representation.",
            "set stick_radius, 0.15, cryoem_structure",
            (
                "set stick_radius, 0.15, "
                "glideem_reference_pose"
            ),
            "set stick_radius, 0.15, top_mcss_poses",
            "set stick_radius, 0.15, binding_site",
            "",
            "orient top_mcss_poses",
            "zoom top_mcss_poses, 8",
            "",
        ]
    )

    return write_pymol_script(
        commands,
        output_dir /
        "pymol_all_mcss_poses_cryo_glideem_ref.pml",
    )

def export_consensus_compound_pose(
    poseviewer: Path,
    compound_id_value: str,
    compound_pose_index: int,
    reference_ligand: Path,
    output_dir: Path,
    metadata: dict,
    binding_site_cutoff: float = 5.0,
) -> Path:
    """
    Export one representative MCSS pose for a consensus-ranked
    compound.

    compound_pose_index is the ordinal occurrence of this
    compound within the PoseViewer, matching the counting used
    during consensus-score collection. It is not the compound's
    rank relative to other compounds in the dock.
    """
    poseviewer = (
        Path(
            poseviewer
        )
        .expanduser()
        .resolve()
    )

    reference_ligand = (
        Path(
            reference_ligand
        )
        .expanduser()
        .resolve()
    )

    output_dir = (
        Path(
            output_dir
        )
        .expanduser()
        .resolve()
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not poseviewer.is_file():
        raise FileNotFoundError(
            f"MCSS PoseViewer not found: {poseviewer}"
        )

    if not reference_ligand.is_file():
        raise FileNotFoundError(
            "GlideEM reference ligand not found: "
            f"{reference_ligand}"
        )

    receptor = None
    selected_pose = None
    compound_pose_counter = 0

    for st in structure.StructureReader(
        str(
            poseviewer
        )
    ):
        if st.property.get(
            "b_glide_receptor",
            False,
        ):
            receptor = st
            continue

        cid = compound_id(
            st
        )

        if cid is None:
            continue

        if (
            str(
                cid
            )
            != str(
                compound_id_value
            )
        ):
            continue

        compound_pose_counter += 1

        if (
            compound_pose_counter
            == compound_pose_index
        ):
            selected_pose = st
            break

    if receptor is None:
        raise RuntimeError(
            f"No receptor CT found in {poseviewer}"
        )

    if selected_pose is None:
        raise RuntimeError(
            "Selected MCSS compound pose was not found: "
            f"compound={compound_id_value}, "
            f"compound_pose_rank={compound_pose_index}, "
            f"PoseViewer={poseviewer}"
        )

    write_structure(
        receptor,
        output_dir
        / "receptor.pdb",
    )

    write_structure(
        selected_pose,
        output_dir
        / "selected_compound_pose.pdb",
    )

    reference_structure = (
        read_first_structure(
            reference_ligand
        )
    )

    write_structure(
        reference_structure,
        output_dir
        / "glideem_reference_pose.pdb",
    )

    metadata_file = (
        output_dir
        / "consensus_pose_metadata.json"
    )

    metadata_file.write_text(
        json.dumps(
            {
                "compound_id": str(
                    compound_id_value
                ),
                "selected_mcss_pose_rank": int(
                    compound_pose_index
                ),
                "poseviewer": str(
                    poseviewer
                ),
                **metadata,
            },
            indent=2,
        )
        + "\n"
    )

    compound_object = (
        "consensus_"
        + safe_pymol_name(
            compound_id_value,
            fallback="compound",
        )
    )

    commands = [
        "reinitialize",
        "bg_color white",
        "",
        "load receptor.pdb, receptor",
        "hide everything",
        "show cartoon, receptor",
        "color gray80, receptor",
        "",
        (
            "load glideem_reference_pose.pdb, "
            "glideem_reference_pose"
        ),
        "show sticks, glideem_reference_pose",
        *element_colour_commands(
            "glideem_reference_pose",
            "salmon",
        ),
        "",
        (
            "load selected_compound_pose.pdb, "
            f"{compound_object}"
        ),
        f"show sticks, {compound_object}",
        *element_colour_commands(
            compound_object,
            "cyan",
        ),
        "",
        (
            "select binding_site, "
            f"byres (receptor within "
            f"{binding_site_cutoff:.1f} of "
            f"({compound_object} or "
            "glideem_reference_pose))"
        ),
        "show sticks, binding_site",
        *element_colour_commands(
            "binding_site",
            "gray80",
        ),
        "",
        "# Consistent thin-stick representation.",
        (
            "set stick_radius, 0.15, "
            "glideem_reference_pose"
        ),
        (
            f"set stick_radius, 0.15, "
            f"{compound_object}"
        ),
        "set stick_radius, 0.15, binding_site",
        "",
        f"orient {compound_object}",
        f"zoom {compound_object}, 8",
        "",
    ]

    return write_pymol_script(
        commands,
        output_dir
        / "view_consensus_pose.pml",
    )

def read_campaign_alignment_transforms(
    alignment_transforms_file: Path,
) -> dict:
    """
    Read campaign receptor-alignment transforms keyed by source name.
    """
    alignment_transforms_file = (
        Path(
            alignment_transforms_file
        )
        .expanduser()
        .resolve()
    )

    if not alignment_transforms_file.is_file():
        raise FileNotFoundError(
            "Campaign receptor-alignment transform table "
            f"not found: {alignment_transforms_file}"
        )

    table = pd.read_csv(
        alignment_transforms_file
    )

    required_columns = {
        "source_name",
        "alignment_reference",
        "rotation_00",
        "rotation_01",
        "rotation_02",
        "rotation_10",
        "rotation_11",
        "rotation_12",
        "rotation_20",
        "rotation_21",
        "rotation_22",
        "translation_x",
        "translation_y",
        "translation_z",
    }

    missing_columns = (
        required_columns
        - set(
            table.columns
        )
    )

    if missing_columns:
        raise RuntimeError(
            "Campaign alignment-transform table is missing "
            f"columns: {sorted(missing_columns)}"
        )

    if table[
        "source_name"
    ].duplicated().any():
        raise RuntimeError(
            "Campaign alignment-transform table contains "
            "duplicate source_name values."
        )

    transforms = {}

    for row in table.itertuples(
        index=False
    ):
        rotation = np.asarray(
            [
                [
                    row.rotation_00,
                    row.rotation_01,
                    row.rotation_02,
                ],
                [
                    row.rotation_10,
                    row.rotation_11,
                    row.rotation_12,
                ],
                [
                    row.rotation_20,
                    row.rotation_21,
                    row.rotation_22,
                ],
            ],
            dtype=float,
        )

        translation = np.asarray(
            [
                row.translation_x,
                row.translation_y,
                row.translation_z,
            ],
            dtype=float,
        )

        if not (
            np.all(
                np.isfinite(
                    rotation
                )
            )
            and np.all(
                np.isfinite(
                    translation
                )
            )
        ):
            raise RuntimeError(
                "Non-finite campaign transform for source "
                f"{row.source_name!r}"
            )

        transforms[
            str(
                row.source_name
            )
        ] = {
            "alignment_reference": str(
                row.alignment_reference
            ),
            "rotation": rotation,
            "translation": translation,
        }

    if not transforms:
        raise RuntimeError(
            "Campaign alignment-transform table is empty."
        )

    return transforms

def export_combined_consensus_poses(
    selected_poses,
    alignment_transforms_file: Path,
    output_dir: Path,
    binding_site_cutoff: float = 5.0,
) -> Path:
    """
    Export one combined PyMOL bundle containing the selected
    consensus pose and corresponding receptor for every compound.

    No GlideEM/cryo-EM reference ligand is included.

    selected_poses must contain:
        consensus_rank
        compound_id
        compound_pose_index
        poseviewer_path
        source_name
        dock_name
    """
    
    output_dir = (
        Path(
            output_dir
        )
        .expanduser()
        .resolve()
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    transforms_by_source = (
        read_campaign_alignment_transforms(
            alignment_transforms_file
        )
    )

    alignment_references = {
        transform[
            "alignment_reference"
        ]
        for transform in (
            transforms_by_source.values()
        )
    }

    if len(
        alignment_references
    ) != 1:
        raise RuntimeError(
            "Consensus export requires all campaign transforms "
            "to use one alignment reference; found: "
            f"{sorted(alignment_references)}"
        )

    alignment_reference = next(
        iter(
            alignment_references
        )
    )    

    commands = [
        "reinitialize",
        "bg_color white",
        "",
    ]

    group_objects = []
    receptor_objects = []
    ligand_objects = []
    export_rows = []
    
    exported_count = 0

    carbon_colours = [
        "cyan",
        "tv_green",
        "violet",
        "orange",
        "marine",
        "magenta",
        "yellow",
        "salmon",
    ]

    for row in selected_poses.itertuples(
        index=False
    ):
        consensus_rank = int(
            row.consensus_rank
        )

        compound_value = str(
            row.compound_id
        )

        safe_id = safe_pymol_name(
            compound_value,
            fallback=(
                f"compound_{consensus_rank:03d}"
            ),
        )

        poseviewer = (
            Path(
                row.poseviewer_path
            )
            .expanduser()
            .resolve()
        )

        if not poseviewer.is_file():
            continue

        compound_pose_index = int(
            row.compound_pose_index
        )

        receptor = None
        selected_pose = None
        compound_pose_counter = 0

        for st in structure.StructureReader(
            str(
                poseviewer
            )
        ):
            if st.property.get(
                "b_glide_receptor",
                False,
            ):
                receptor = st
                continue

            cid = compound_id(
                st
            )

            if cid is None:
                continue

            if (
                str(
                    cid
                )
                != compound_value
            ):
                continue

            compound_pose_counter += 1

            if (
                compound_pose_counter
                == compound_pose_index
            ):
                selected_pose = st
                break

                # ----------------------------------------------
        # Retrieve this structure's campaign transform
        # ----------------------------------------------

        source_name = str(
            row.source_name
        )

        if (
            source_name
            not in transforms_by_source
        ):
            raise RuntimeError(
                "No campaign alignment transform was found for "
                f"source structure {source_name!r}. "
                "The consensus ranking and campaign transform "
                "table may come from different runs."
            )

        source_transform = (
            transforms_by_source[
                source_name
            ]
        )

        rotation = source_transform[
            "rotation"
        ]

        translation = source_transform[
            "translation"
        ]

        # ----------------------------------------------
        # Apply exactly the same rigid transform to the
        # MCSS receptor and its selected ligand.
        # ----------------------------------------------

        aligned_receptor = (
            transform_structure_copy(
                input_structure=receptor,
                rotation=rotation,
                translation=translation,
            )
        )

        aligned_selected_pose = (
            transform_structure_copy(
                input_structure=selected_pose,
                rotation=rotation,
                translation=translation,
            )
        )

        aligned_receptor.title = (
            f"{source_name}_aligned_to_"
            f"{alignment_reference}"
        )

        aligned_selected_pose.title = (
            f"consensus_rank{consensus_rank:03d}_"
            f"{safe_id}_aligned_to_"
            f"{alignment_reference}"
        )

        # ----------------------------------------------
        # Construct unique filenames / PyMOL object names
        # after the source structure is known.
        # ----------------------------------------------

        safe_source = safe_pymol_name(
            source_name,
            fallback="source",
        )

        prefix = (
            f"rank{consensus_rank:03d}_"
            f"{safe_id}_"
            f"{safe_source}"
        )

        receptor_object = (
            f"receptor_{prefix}"
        )

        ligand_object = (
            f"ligand_{prefix}"
        )

        site_object = (
            f"site_{prefix}"
        )

        group_object = (
            f"consensus_{prefix}"
        )

        group_objects.append(
            group_object
        )

        receptor_objects.append(
            receptor_object
        )

        ligand_objects.append(
            ligand_object
        )

        receptor_file = (
            output_dir
            / f"{receptor_object}.pdb"
        )

        ligand_file = (
            output_dir
            / f"{ligand_object}.pdb"
        )

        # ----------------------------------------------
        # Write the transformed structures, not the native
        # PoseViewer coordinates.
        # ----------------------------------------------

        write_structure(
            aligned_receptor,
            receptor_file,
        )

        write_structure(
            aligned_selected_pose,
            ligand_file,
        )        

        carbon_colour = (
            carbon_colours[
                exported_count
                % len(
                    carbon_colours
                )
            ]
        )

        commands.extend(
            [
                (
                    f"load {receptor_file.name}, "
                    f"{receptor_object}"
                ),
                (
                    f"hide everything, "
                    f"{receptor_object}"
                ),
                (
                    f"show cartoon, "
                    f"{receptor_object}"
                ),
                (
                    f"color gray80, "
                    f"{receptor_object}"
                ),
                "",
                (
                    f"load {ligand_file.name}, "
                    f"{ligand_object}"
                ),
                (
                    f"show sticks, "
                    f"{ligand_object}"
                ),
                *element_colour_commands(
                    ligand_object,
                    carbon_colour,
                ),
                "",
                (
                    f"select {site_object}, "
                    f"byres ({receptor_object} "
                    f"within "
                    f"{binding_site_cutoff:.1f} "
                    f"of {ligand_object})"
                ),
                (
                    f"show sticks, "
                    f"{site_object}"
                ),
                *element_colour_commands(
                    site_object,
                    "gray80",
                ),
                "",
                (
                    f"set stick_radius, 0.15, "
                    f"{ligand_object}"
                ),
                (
                    f"set stick_radius, 0.15, "
                    f"{site_object}"
                ),
                "",
                (
                    f"group {group_object}, "
                    f"{receptor_object} "
                    f"{ligand_object} "
                    f"{site_object}"
                ),
                "",
            ]
        )

        export_rows.append(
            {
                "consensus_rank": (
                    consensus_rank
                ),
                "compound_id": (
                    compound_value
                ),
                "source_name": (
                    source_name
                ),
                "dock_name": str(
                    row.dock_name
                ),
                "compound_pose_index": (
                    compound_pose_index
                ),
                "alignment_reference": (
                    alignment_reference
                ),
                "receptor_object": (
                    receptor_object
                ),
                "ligand_object": (
                    ligand_object
                ),
                "receptor_file": (
                    receptor_file.name
                ),
                "ligand_file": (
                    ligand_file.name
                ),
            }
        )        

        exported_count += 1

    if exported_count == 0:
        raise RuntimeError(
            "No consensus receptor/ligand pairs could "
            "be exported."
        )

    first_group = group_objects[
        0
    ]

    first_ligand = ligand_objects[
        0
    ]

    commands.extend(
        [
            "",
            "# Start with the highest-ranked compound visible.",
        ]
    )

    commands.extend(
        [
            "",
            (
                "# All selected compounds are shown in the common "
                "campaign coordinate frame."
            ),
            (
                "# Only the receptor associated with consensus "
                "rank 1 is shown initially."
            ),
        ]
    )

    for receptor_object in receptor_objects[
        1:
    ]:
        commands.append(
            f"disable {receptor_object}"
        )

    commands.extend(
        [
            "",
            (
                "# All ranked receptor/pose pairs are retained "
                "as individual PyMOL groups."
            ),
            "",
            (
                "select consensus_ligands, "
                + " or ".join(
                    ligand_objects
                )
            ),
            "orient consensus_ligands",
            "zoom consensus_ligands, 8",
            "",
        ]
    )

    pd.DataFrame(
        export_rows
    ).to_csv(
        output_dir
        / "combined_consensus_pose_manifest.csv",
        index=False,
    )    

    script_file = write_pymol_script(
        commands,
        output_dir
        / "view_all_consensus_poses.pml",
    )

    return script_file