from __future__ import annotations

import csv
import re

from collections import defaultdict
from pathlib import Path

from schrodinger import structure


ELEMENT_COLOURS = {
    "N": "blue",
    "O": "red",
    "S": "yellow",
    "P": "orange",
    "H": "white",
    "F": "green",
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


def element_colour_commands(
    object_name: str,
    carbon_colour: str,
) -> list[str]:
    """
    Return PyMOL commands that colour carbon atoms using
    the object-specific colour and other atoms by element.
    """

    commands = [
        (
            f"color {carbon_colour}, "
            f"({object_name} and elem C)"
        ),
    ]

    for element, colour in (
        ELEMENT_COLOURS.items()
    ):
        commands.append(
            f"color {colour}, "
            f"({object_name} and elem {element})"
        )

    return commands


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


def export_glideem_pose_ensemble(
    pose_structures,
    cluster_by_pose,
    cryoem_reference: Path,
    output_dir: Path,
    job_name: str,
) -> Path:
    output_dir = (
        Path(output_dir)
        .expanduser()
        .resolve()
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    cryoem_structure = read_first_structure(
        cryoem_reference
    )

    write_structure(
        cryoem_structure,
        output_dir /
        "cryoem_structure.pdb",
    )

    exported = []

    for record in pose_structures:
        pose_number = int(
            record["pose_number"]
        )

        cluster_id = int(
            cluster_by_pose[
                pose_number
            ]
        )

        object_name = (
            f"pose{pose_number:02d}"
            f"_cluster{cluster_id:02d}"
        )

        pose_file = write_structure(
            record["structure"],
            output_dir /
            f"{object_name}.pdb",
        )

        exported.append(
            {
                "pose_number": pose_number,
                "cluster_id": cluster_id,
                "object_name": object_name,
                "pose_file": pose_file,
            }
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

    commands = [
        "reinitialize",
        "bg_color white",
        "",
        "load cryoem_structure.pdb, cryoem_structure",
        "show sticks, cryoem_structure",
        *element_colour_commands(
            "cryoem_structure",
            "gray50",
        ),
        "",
    ]

    cluster_objects = defaultdict(
        list
    )

    for record in exported:
        cluster_id = record[
            "cluster_id"
        ]

        object_name = record[
            "object_name"
        ]

        pose_file = record[
            "pose_file"
        ]

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

    for cluster_name, objects in sorted(
        cluster_objects.items()
    ):
        commands.append(
            f"group {cluster_name}, "
            + " ".join(objects)
        )

    commands.extend(
        [
            "",
            "group glideem_pose_clusters, "
            + " ".join(
                sorted(
                    cluster_objects
                )
            ),
            "",
            "orient glideem_pose_clusters",
            "zoom glideem_pose_clusters, 6",
            "",
        ]
    )

    return write_pymol_script(
        commands,
        output_dir /
        "load_glideem_ensemble.pml",
    )


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

    if ligand_objects:
        commands.append(
            "group top_mcss_poses, "
            + " ".join(
                ligand_objects
            )
        )

    ligand_selection = (
        "cryoem_structure or "
        "glideem_reference_pose or "
        "top_mcss_poses"
    )

    commands.extend(
        [
            "",
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