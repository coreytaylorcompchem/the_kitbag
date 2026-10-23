from __future__ import annotations

import re
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from scipy.cluster.hierarchy import (
    dendrogram,
    fcluster,
    linkage,
)
from scipy.ndimage import gaussian_filter
from scipy.spatial.distance import (
    cdist,
    pdist,
    jensenshannon,
    squareform,
)

from schrodinger import structure

try:
    import mrcfile

except ImportError:
    mrcfile = None

from modules.utils.glideem_helpers import (
    find_glideem_ligand_indices,
    ligand_heavy_atom_coordinates,
)

from modules.utils.constants import THREE_TO_ONE, PYMOL_COLOURS, _safe_object_name

# ============================================================
# Basic numerical utilities
# ============================================================

def _normalise_probability(values: np.ndarray) -> np.ndarray:
    """
    Return a finite probability vector.

    A zero vector remains a zero vector rather than producing NaNs.
    """
    values = np.asarray(
        values,
        dtype=np.float64,
    )

    values = np.where(
        np.isfinite(values),
        values,
        0.0,
    )

    values = np.clip(
        values,
        0.0,
        None,
    )

    total = float(
        values.sum()
    )

    if total <= 0.0:
        return np.zeros_like(
            values,
            dtype=np.float64,
        )

    return values / total


def _shannon_entropy(
    probabilities: np.ndarray,
) -> float:
    """
    Unnormalised Shannon entropy.

    exp(entropy) gives the effective number of occupied families.
    """
    probabilities = np.asarray(
        probabilities,
        dtype=float,
    )

    probabilities = probabilities[
        probabilities > 0.0
    ]

    if probabilities.size == 0:
        return 0.0

    return float(
        -np.sum(
            probabilities
            * np.log(probabilities)
        )
    )


def _normalised_entropy(
    probabilities: np.ndarray,
) -> float:
    probabilities = np.asarray(
        probabilities,
        dtype=float,
    )

    probabilities = probabilities[
        probabilities > 0.0
    ]

    if probabilities.size <= 1:
        return 0.0

    entropy = _shannon_entropy(
        probabilities
    )

    return float(
        entropy
        / math.log(
            probabilities.size
        )
    )


# ============================================================
# Protein alignment
# ============================================================

def _atom_name(atom) -> str:
    return (
        getattr(
            atom,
            "pdbname",
            "",
        )
        or ""
    ).strip().upper()


def _residue_key(atom) -> tuple:
    """
    Build an alignment key from PDB identifiers.

    This assumes equivalent receptor residues retain compatible
    chain/residue numbering across input structures.
    """
    chain = (
        getattr(
            atom,
            "chain",
            "",
        )
        or ""
    ).strip()

    resnum = int(
        getattr(
            atom,
            "resnum",
            0,
        )
    )

    inscode = (
        getattr(
            atom,
            "inscode",
            "",
        )
        or ""
    ).strip()

    pdbres = (
        getattr(
            atom,
            "pdbres",
            "",
        )
        or ""
    ).strip().upper()

    return (
        chain,
        resnum,
        inscode,
        pdbres,
    )


def _all_ca_coordinates(
    receptor_st,
    atom_name: str = "CA",
) -> dict[tuple, np.ndarray]:
    """
    Return receptor alignment atoms keyed by PDB residue identity.
    """
    coordinates = {}

    for atom in receptor_st.atom:
        if _atom_name(atom) != atom_name.upper():
            continue

        coordinates[
            _residue_key(atom)
        ] = np.array(
            [
                atom.x,
                atom.y,
                atom.z,
            ],
            dtype=float,
        )

    return coordinates


def _minimum_distance_to_points(
    points: np.ndarray,
    reference_points: np.ndarray,
) -> np.ndarray:
    if not len(points):
        return np.array(
            [],
            dtype=float,
        )

    if not len(reference_points):
        return np.full(
            len(points),
            np.inf,
            dtype=float,
        )

    distances = cdist(
        points,
        reference_points,
    )

    return distances.min(
        axis=1
    )


def pocket_alignment_coordinates(
    receptor_st,
    reference_ligand_coordinates: np.ndarray,
    pocket_cutoff: float,
    atom_name: str = "CA",
) -> dict[tuple, np.ndarray]:
    """
    Select receptor alignment atoms whose coordinates lie within
    pocket_cutoff of the experimental reference ligand.

    Selection is performed separately for every structure. The final
    alignment uses the residue-key intersection between structures.
    """
    all_coordinates = _all_ca_coordinates(
        receptor_st,
        atom_name=atom_name,
    )

    if not all_coordinates:
        return {}

    keys = list(
        all_coordinates
    )

    coordinates = np.asarray(
        [
            all_coordinates[key]
            for key in keys
        ],
        dtype=float,
    )

    minimum_distances = (
        _minimum_distance_to_points(
            coordinates,
            reference_ligand_coordinates,
        )
    )

    return {
        key: coordinate
        for key, coordinate, distance in zip(
            keys,
            coordinates,
            minimum_distances,
        )
        if distance <= pocket_cutoff
    }


def kabsch_transform(
    mobile_coordinates: np.ndarray,
    reference_coordinates: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Return rotation and translation that transform mobile row-vector
    coordinates into the reference coordinate frame:

        aligned = mobile @ rotation + translation
    """
    mobile_coordinates = np.asarray(
        mobile_coordinates,
        dtype=float,
    )

    reference_coordinates = np.asarray(
        reference_coordinates,
        dtype=float,
    )

    if (
        mobile_coordinates.ndim != 2
        or reference_coordinates.ndim != 2
        or mobile_coordinates.shape[1] != 3
        or reference_coordinates.shape[1] != 3
    ):
        raise ValueError(
            "Kabsch input coordinates must have shape (n, 3)."
        )

    if (
        mobile_coordinates.shape
        != reference_coordinates.shape
    ):
        raise ValueError(
            "Mobile and reference alignment coordinates "
            "must have the same shape."
        )

    mobile_centroid = (
        mobile_coordinates.mean(
            axis=0
        )
    )

    reference_centroid = (
        reference_coordinates.mean(
            axis=0
        )
    )

    mobile_centered = (
        mobile_coordinates
        - mobile_centroid
    )

    reference_centered = (
        reference_coordinates
        - reference_centroid
    )

    covariance = (
        mobile_centered.T
        @ reference_centered
    )

    u_matrix, _, vt_matrix = (
        np.linalg.svd(
            covariance
        )
    )

    rotation = (
        u_matrix
        @ vt_matrix
    )

    # Prevent reflection.
    if np.linalg.det(rotation) < 0:
        u_matrix[:, -1] *= -1.0
        rotation = (
            u_matrix
            @ vt_matrix
        )

    translation = (
        reference_centroid
        - mobile_centroid @ rotation
    )

    return (
        rotation,
        translation,
    )


def apply_transform(
    coordinates: np.ndarray,
    rotation: np.ndarray,
    translation: np.ndarray,
) -> np.ndarray:
    coordinates = np.asarray(
        coordinates,
        dtype=float,
    )

    return (
        coordinates
        @ rotation
        + translation
    )


def alignment_rmsd(
    mobile_coordinates: np.ndarray,
    reference_coordinates: np.ndarray,
    rotation: np.ndarray,
    translation: np.ndarray,
) -> float:
    aligned = apply_transform(
        mobile_coordinates,
        rotation,
        translation,
    )

    squared_distances = np.sum(
        (
            aligned
            - reference_coordinates
        ) ** 2,
        axis=1,
    )

    return float(
        np.sqrt(
            np.mean(
                squared_distances
            )
        )
    )


# ============================================================
# Structure and pose loading
# ============================================================

def read_first_structure(
    input_file: Path,
):
    entries = list(
        structure.StructureReader(
            str(input_file)
        )
    )

    if not entries:
        raise RuntimeError(
            f"No structures found in {input_file}"
        )

    return entries[0]


def read_reference_ligand_coordinates(
    reference_file: Path,
) -> np.ndarray:
    ligand_st = read_first_structure(
        reference_file
    )

    coordinates = []

    for atom in ligand_st.atom:
        if atom.element.upper() == "H":
            continue

        coordinates.append(
            [
                atom.x,
                atom.y,
                atom.z,
            ]
        )

    if not coordinates:
        raise RuntimeError(
            f"No reference ligand heavy atoms found in "
            f"{reference_file}"
        )

    return np.asarray(
        coordinates,
        dtype=float,
    )

def maximum_internal_distance_error(
    original_coordinates: np.ndarray,
    transformed_coordinates: np.ndarray,
) -> float:
    original_distances = cdist(
        original_coordinates,
        original_coordinates,
    )

    transformed_distances = cdist(
        transformed_coordinates,
        transformed_coordinates,
    )

    return float(
        np.max(
            np.abs(
                original_distances
                - transformed_distances
            )
        )
    )


def load_aligned_glideem_poses(
    archive: Path,
    ligand_resname: str,
    ligand_from_smiles: bool,
    rotation: np.ndarray,
    translation: np.ndarray,
) -> list[dict]:
    """
    Load every usable GlideEM ligand pose and transform it into the
    campaign reference frame.

    Heavy-atom coordinates are retained for numerical analysis.
    All ligand atoms, including hydrogens, are transformed for export.
    """
    poses = []

    for pose_number, st in enumerate(
        structure.StructureReader(
            str(archive)
        ),
        start=1,
    ):
        ligand_indices = (
            find_glideem_ligand_indices(
                st,
                ligand_resname=ligand_resname,
                ligand_from_smiles=ligand_from_smiles,
            )
        )

        if not ligand_indices:
            continue

        # Numerical analysis uses the existing heavy-atom helper.
        heavy_coordinates, elements = (
            ligand_heavy_atom_coordinates(
                st,
                ligand_indices,
            )
        )

        aligned_heavy_coordinates = (
            apply_transform(
                heavy_coordinates,
                rotation,
                translation,
            )
        )

        # Export structure contains the complete ligand.
        ligand_st = st.extract(
            ligand_indices
        )

        # Transform every atom in the extracted ligand in its own
        # original position. Do not assign heavy-atom coordinates
        # sequentially to a structure that may contain hydrogens.
        all_atom_coordinates = np.asarray(
            [
                [
                    atom.x,
                    atom.y,
                    atom.z,
                ]
                for atom in ligand_st.atom
            ],
            dtype=float,
        )

        aligned_all_atom_coordinates = (
            apply_transform(
                all_atom_coordinates,
                rotation,
                translation,
            )
        )

        for atom, coordinate in zip(
            ligand_st.atom,
            aligned_all_atom_coordinates,
        ):
            atom.x = float(
                coordinate[0]
            )
            atom.y = float(
                coordinate[1]
            )
            atom.z = float(
                coordinate[2]
            )

        distance_error = (
            maximum_internal_distance_error(
                all_atom_coordinates,
                aligned_all_atom_coordinates,
            )
        )

        if distance_error > 1.0e-5:
            raise RuntimeError(
                f"Rigid transformation changed internal distances "
                f"for pose {pose_number}: "
                f"maximum error={distance_error:.6g} A"
            )

        poses.append(
            {
                "pose_number": pose_number,
                "coordinates": (
                    aligned_heavy_coordinates
                ),
                "elements": elements,
                "structure": ligand_st,
                "centroid": (
                    aligned_heavy_coordinates.mean(
                        axis=0
                    )
                ),
            }
        )

    return poses


# ============================================================
# Common 3D grid and densities
# ============================================================

def make_common_grid(
    all_coordinates: list[np.ndarray],
    spacing: float,
    padding: float,
) -> dict:
    """
    Construct one common orthogonal grid containing every aligned pose.
    """
    if spacing <= 0.0:
        raise ValueError(
            "Density grid spacing must be positive."
        )

    if not all_coordinates:
        raise ValueError(
            "Cannot build a common grid without coordinates."
        )

    stacked = np.vstack(
        all_coordinates
    )

    minimum = (
        stacked.min(
            axis=0
        )
        - padding
    )

    maximum = (
        stacked.max(
            axis=0
        )
        + padding
    )

    shape = (
        np.ceil(
            (
                maximum
                - minimum
            )
            / spacing
        ).astype(int)
        + 1
    )

    if np.any(shape < 2):
        raise RuntimeError(
            f"Invalid common density grid shape: {shape}"
        )

    return {
        "origin": minimum.astype(float),
        "spacing": float(spacing),
        "shape": tuple(
            int(value)
            for value in shape
        ),
        "maximum": maximum.astype(float),
    }


def coordinates_to_grid_indices(
    coordinates: np.ndarray,
    grid: dict,
) -> np.ndarray:
    origin = np.asarray(
        grid["origin"],
        dtype=float,
    )

    spacing = float(
        grid["spacing"]
    )

    indices = np.rint(
        (
            coordinates
            - origin
        )
        / spacing
    ).astype(int)

    shape = np.asarray(
        grid["shape"],
        dtype=int,
    )

    valid = (
        np.all(
            indices >= 0,
            axis=1,
        )
        & np.all(
            indices < shape,
            axis=1,
        )
    )

    return indices[
        valid
    ]


def pose_probability_density(
    coordinates: np.ndarray,
    grid: dict,
    gaussian_sigma_angstrom: float,
) -> np.ndarray:
    """
    Convert one pose to a normalized heavy-atom probability density.

    Normalizing every pose to unit mass prevents larger ligands from
    automatically dominating pooled structure and target maps.
    """
    density = np.zeros(
        grid["shape"],
        dtype=np.float32,
    )

    indices = coordinates_to_grid_indices(
        coordinates,
        grid,
    )

    if not len(indices):
        return density

    for x_index, y_index, z_index in indices:
        density[
            x_index,
            y_index,
            z_index,
        ] += 1.0

    sigma_voxels = (
        gaussian_sigma_angstrom
        / float(
            grid["spacing"]
        )
    )

    if sigma_voxels > 0.0:
        density = gaussian_filter(
            density,
            sigma=sigma_voxels,
            mode="constant",
        ).astype(
            np.float32
        )

    total = float(
        density.sum()
    )

    if total > 0.0:
        density /= total

    return density


def pooled_probability_density(
    densities: list[np.ndarray],
) -> np.ndarray:
    if not densities:
        raise ValueError(
            "Cannot pool an empty density collection."
        )

    pooled = np.mean(
        np.stack(
            densities,
            axis=0,
        ),
        axis=0,
    )

    total = float(
        pooled.sum()
    )

    if total > 0.0:
        pooled /= total

    return pooled.astype(
        np.float32
    )


def density_variability(
    densities: list[np.ndarray],
) -> np.ndarray:
    if not densities:
        raise ValueError(
            "Cannot calculate variability for an empty collection."
        )

    return np.std(
        np.stack(
            densities,
            axis=0,
        ),
        axis=0,
    ).astype(
        np.float32
    )


def density_difference(
    first: np.ndarray,
    second: np.ndarray,
) -> np.ndarray:
    return (
        np.asarray(
            first,
            dtype=np.float32,
        )
        - np.asarray(
            second,
            dtype=np.float32,
        )
    )


# ============================================================
# Global pose-family clustering
# ============================================================

def pairwise_density_distance_matrix(
    densities: list[np.ndarray],
    metric: str = "jensenshannon",
) -> np.ndarray:
    """
    Compare differently sized ligands by their normalized 3D density
    distributions rather than atom-to-atom RMSD.
    """
    n_items = len(
        densities
    )

    matrix = np.zeros(
        (
            n_items,
            n_items,
        ),
        dtype=np.float32,
    )

    flattened = [
        _normalise_probability(
            density.ravel()
        )
        for density in densities
    ]

    for first_index in range(
        n_items
    ):
        for second_index in range(
            first_index + 1,
            n_items,
        ):
            first = flattened[
                first_index
            ]

            second = flattened[
                second_index
            ]

            if metric == "jensenshannon":
                distance = float(
                    jensenshannon(
                        first,
                        second,
                        base=2.0,
                    )
                )

            else:
                raise ValueError(
                    "Unsupported global pose-density "
                    f"distance metric: {metric}"
                )

            if not np.isfinite(
                distance
            ):
                distance = 1.0

            matrix[
                first_index,
                second_index,
            ] = distance

            matrix[
                second_index,
                first_index,
            ] = distance

    return matrix


def cluster_global_pose_families(
    distance_matrix: np.ndarray,
    distance_cutoff: float,
    linkage_method: str = "average",
) -> tuple[np.ndarray, np.ndarray | None]:
    n_items = int(
        distance_matrix.shape[0]
    )

    if n_items == 0:
        raise ValueError(
            "No poses supplied for global clustering."
        )

    if n_items == 1:
        return (
            np.ones(
                1,
                dtype=int,
            ),
            None,
        )

    condensed = squareform(
        distance_matrix,
        checks=False,
    )

    linkage_matrix = linkage(
        condensed,
        method=linkage_method,
    )

    family_ids = fcluster(
        linkage_matrix,
        t=distance_cutoff,
        criterion="distance",
    )

    # Re-number by population so Family 1 is the largest family.
    unique, counts = np.unique(
        family_ids,
        return_counts=True,
    )

    population_order = [
        int(value)
        for value in unique[
            np.argsort(
                -counts
            )
        ]
    ]

    remapping = {
        old_id: new_id
        for new_id, old_id in enumerate(
            population_order,
            start=1,
        )
    }

    family_ids = np.asarray(
        [
            remapping[
                int(family_id)
            ]
            for family_id in family_ids
        ],
        dtype=int,
    )

    return (
        family_ids,
        linkage_matrix,
    )


def family_medoid_indices(
    distance_matrix: np.ndarray,
    family_ids: np.ndarray,
) -> dict[int, int]:
    """
    Return the global pose index minimizing mean within-family distance.
    """
    medoids = {}

    for family_id in sorted(
        set(
            family_ids.tolist()
        )
    ):
        member_indices = np.where(
            family_ids == family_id
        )[0]

        if len(member_indices) == 1:
            medoids[
                int(family_id)
            ] = int(
                member_indices[0]
            )
            continue

        submatrix = distance_matrix[
            np.ix_(
                member_indices,
                member_indices,
            )
        ]

        mean_distances = (
            submatrix.mean(
                axis=1
            )
        )

        local_medoid = int(
            np.argmin(
                mean_distances
            )
        )

        medoids[
            int(family_id)
        ] = int(
            member_indices[
                local_medoid
            ]
        )

    return medoids


# ============================================================
# Distributional summaries
# ============================================================

def structure_family_composition(
    pose_table: pd.DataFrame,
) -> pd.DataFrame:
    counts = (
        pose_table
        .groupby(
            [
                "source_name",
                "global_family_id",
            ],
            dropna=False,
        )
        .size()
        .rename(
            "n_poses"
        )
        .reset_index()
    )

    totals = (
        counts
        .groupby(
            "source_name"
        )[
            "n_poses"
        ]
        .transform(
            "sum"
        )
    )

    counts[
        "family_fraction"
    ] = (
        counts[
            "n_poses"
        ]
        / totals
    )

    return counts


def structure_ambiguity_summary(
    composition: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    for source_name, group in composition.groupby(
        "source_name"
    ):
        probabilities = (
            group[
                "family_fraction"
            ]
            .to_numpy(
                dtype=float
            )
        )

        entropy = _shannon_entropy(
            probabilities
        )

        rows.append(
            {
                "source_name": source_name,
                "n_global_families": int(
                    np.sum(
                        probabilities > 0.0
                    )
                ),
                "dominant_family_fraction": float(
                    probabilities.max()
                ),
                "family_entropy": entropy,
                "normalised_family_entropy": (
                    _normalised_entropy(
                        probabilities
                    )
                ),
                "effective_n_families": float(
                    np.exp(
                        entropy
                    )
                ),
            }
        )

    return pd.DataFrame(
        rows
    )


def family_support_summary(
    pose_table: pd.DataFrame,
) -> pd.DataFrame:
    total_structures = int(
        pose_table[
            "source_name"
        ].nunique()
    )

    rows = []

    for family_id, group in pose_table.groupby(
        "global_family_id"
    ):
        n_structures = int(
            group[
                "source_name"
            ].nunique()
        )

        rows.append(
            {
                "global_family_id": int(
                    family_id
                ),
                "n_poses": int(
                    len(
                        group
                    )
                ),
                "n_supporting_structures": n_structures,
                "structure_support_fraction": float(
                    n_structures
                    / total_structures
                    if total_structures
                    else np.nan
                ),
            }
        )

    return (
        pd.DataFrame(
            rows
        )
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
        .reset_index(
            drop=True
        )
    )


def structure_js_distance_matrix(
    structure_densities: dict[str, np.ndarray],
) -> pd.DataFrame:
    names = sorted(
        structure_densities
    )

    matrix = np.zeros(
        (
            len(names),
            len(names),
        ),
        dtype=float,
    )

    flattened = {
        name: _normalise_probability(
            structure_densities[
                name
            ].ravel()
        )
        for name in names
    }

    for first_index, first_name in enumerate(
        names
    ):
        for second_index in range(
            first_index + 1,
            len(names),
        ):
            second_name = names[
                second_index
            ]

            distance = float(
                jensenshannon(
                    flattened[
                        first_name
                    ],
                    flattened[
                        second_name
                    ],
                    base=2.0,
                )
            )

            matrix[
                first_index,
                second_index,
            ] = distance

            matrix[
                second_index,
                first_index,
            ] = distance

    return pd.DataFrame(
        matrix,
        index=names,
        columns=names,
    )

def transform_structure_copy(
    input_structure,
    rotation: np.ndarray,
    translation: np.ndarray,
):
    transformed = input_structure.copy()

    coordinates = np.asarray(
        [
            [
                atom.x,
                atom.y,
                atom.z,
            ]
            for atom in transformed.atom
        ],
        dtype=float,
    )

    transformed_coordinates = apply_transform(
        coordinates,
        rotation,
        translation,
    )

    for atom, coordinate in zip(
        transformed.atom,
        transformed_coordinates,
    ):
        atom.x = float(
            coordinate[0]
        )
        atom.y = float(
            coordinate[1]
        )
        atom.z = float(
            coordinate[2]
        )

    return transformed

def centroid_dispersion_summary(
    pose_table: pd.DataFrame,
) -> pd.DataFrame:
    """
    Decompose simple centroid variance into within-structure and
    between-structure components.

    This is descriptive rather than a fitted mixed-effects model.
    """
    coordinate_columns = [
        "centroid_x",
        "centroid_y",
        "centroid_z",
    ]

    all_centroids = (
        pose_table[
            coordinate_columns
        ]
        .to_numpy(
            dtype=float
        )
    )

    global_centroid = (
        all_centroids.mean(
            axis=0
        )
    )

    total_ms_distance = float(
        np.mean(
            np.sum(
                (
                    all_centroids
                    - global_centroid
                ) ** 2,
                axis=1,
            )
        )
    )

    structure_rows = []
    structure_centroids = []
    structure_weights = []
    pooled_within_sum = 0.0
    pooled_within_n = 0

    for source_name, group in pose_table.groupby(
        "source_name"
    ):
        centroids = (
            group[
                coordinate_columns
            ]
            .to_numpy(
                dtype=float
            )
        )

        centre = (
            centroids.mean(
                axis=0
            )
        )

        squared_distances = np.sum(
            (
                centroids
                - centre
            ) ** 2,
            axis=1,
        )

        mean_squared_dispersion = float(
            np.mean(
                squared_distances
            )
        )

        structure_rows.append(
            {
                "source_name": source_name,
                "n_poses": int(
                    len(
                        centroids
                    )
                ),
                "centroid_x": float(
                    centre[0]
                ),
                "centroid_y": float(
                    centre[1]
                ),
                "centroid_z": float(
                    centre[2]
                ),
                "within_structure_centroid_rms": float(
                    np.sqrt(
                        mean_squared_dispersion
                    )
                ),
            }
        )

        structure_centroids.append(
            centre
        )

        structure_weights.append(
            len(
                centroids
            )
        )

        pooled_within_sum += float(
            squared_distances.sum()
        )

        pooled_within_n += int(
            len(
                centroids
            )
        )

    structure_centroids = np.asarray(
        structure_centroids,
        dtype=float,
    )

    structure_weights = np.asarray(
        structure_weights,
        dtype=float,
    )

    between_ms_distance = float(
        np.average(
            np.sum(
                (
                    structure_centroids
                    - global_centroid
                ) ** 2,
                axis=1,
            ),
            weights=structure_weights,
        )
    )

    within_ms_distance = float(
        pooled_within_sum
        / pooled_within_n
        if pooled_within_n
        else 0.0
    )

    overall = pd.DataFrame(
        [
            {
                "total_centroid_rms": float(
                    np.sqrt(
                        total_ms_distance
                    )
                ),
                "within_structure_centroid_rms": float(
                    np.sqrt(
                        within_ms_distance
                    )
                ),
                "between_structure_centroid_rms": float(
                    np.sqrt(
                        between_ms_distance
                    )
                ),
                "within_variance_fraction": float(
                    within_ms_distance
                    / total_ms_distance
                    if total_ms_distance > 0
                    else np.nan
                ),
                "between_variance_fraction": float(
                    between_ms_distance
                    / total_ms_distance
                    if total_ms_distance > 0
                    else np.nan
                ),
            }
        ]
    )

    return (
        pd.DataFrame(
            structure_rows
        ),
        overall,
    )


# ============================================================
# Plotting
# ============================================================

def save_structure_family_heatmap(
    composition: pd.DataFrame,
    ambiguity: pd.DataFrame,
    output_file: Path,
):
    matrix = (
        composition
        .pivot(
            index="source_name",
            columns="global_family_id",
            values="family_fraction",
        )
        .fillna(
            0.0
        )
    )

    ambiguity_by_structure = (
        ambiguity
        .set_index(
            "source_name"
        )
    )

    matrix = matrix.loc[
        ambiguity.sort_values(
            "effective_n_families",
            ascending=False,
        )[
            "source_name"
        ]
    ]

    figure_width = min(
        20.0,
        max(
            8.0,
            0.7
            * matrix.shape[1]
            + 4.0,
        ),
    )

    figure_height = max(
        4.5,
        0.45
        * matrix.shape[0]
        + 2.0,
    )

    figure, axis = plt.subplots(
        figsize=(
            figure_width,
            figure_height,
        )
    )

    matrix_values = matrix.to_numpy(
        dtype=float
    )

    observed_maximum = float(
        np.nanmax(
            matrix_values
        )
    )

    colour_maximum = max(
        observed_maximum,
        0.10,
    )

    image = axis.imshow(
        matrix_values,
        aspect="auto",
        vmin=0.0,
        vmax=colour_maximum,
        cmap="Blues",
    )

    axis.set_xticks(
        np.arange(
            matrix.shape[1]
        )
    )

    axis.set_xticklabels(
        [
            f"Family {int(value)}"
            for value in matrix.columns
        ],
        rotation=45,
        ha="right",
    )

    labels = []

    for source_name in matrix.index:
        effective_n = float(
            ambiguity_by_structure.loc[
                source_name,
                "effective_n_families",
            ]
        )

        labels.append(
            f"{source_name}  "
            f"(N_eff={effective_n:.2f})"
        )

    axis.set_yticks(
        np.arange(
            matrix.shape[0]
        )
    )

    axis.set_yticklabels(
        labels
    )

    axis.set_xlabel(
        "Global 3D pose family"
    )

    axis.set_ylabel(
        "Cryo-EM structure / density"
    )

    axis.set_title(
        "GlideEM campaign binding-mode landscape"
    )

    colour_bar = figure.colorbar(
        image,
        ax=axis,
    )

    colour_bar.set_label(
        "Fraction of structure's GlideEM poses"
    )

    figure.tight_layout()

    figure.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(
        figure
    )


def save_structure_js_heatmap(
    js_matrix: pd.DataFrame,
    output_file: Path,
):
    figure_size = max(
        6.0,
        0.55
        * len(
            js_matrix
        )
        + 2.0,
    )

    figure, axis = plt.subplots(
        figsize=(
            figure_size,
            figure_size,
        )
    )

    image = axis.imshow(
        js_matrix.to_numpy(
            dtype=float
        ),
        cmap="magma",
        vmin=0.0,
        vmax=max(
            1.0,
            float(
                np.nanmax(
                    js_matrix.to_numpy(
                        dtype=float
                    )
                )
            ),
        ),
    )

    axis.set_xticks(
        np.arange(
            len(
                js_matrix.columns
            )
        )
    )

    axis.set_xticklabels(
        js_matrix.columns,
        rotation=45,
        ha="right",
    )

    axis.set_yticks(
        np.arange(
            len(
                js_matrix.index
            )
        )
    )

    axis.set_yticklabels(
        js_matrix.index
    )

    axis.set_title(
        "Between-structure 3D landscape divergence"
    )

    colour_bar = figure.colorbar(
        image,
        ax=axis,
    )

    colour_bar.set_label(
        "Jensen-Shannon distance"
    )

    figure.tight_layout()

    figure.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(
        figure
    )


def save_global_family_dendrogram(
    linkage_matrix: np.ndarray | None,
    labels: list[str],
    cutoff: float,
    output_file: Path,
    distance_label: str = "Distance",
    title: str = "Global GlideEM pose-family clustering",
):
    if linkage_matrix is None:
        return

    figure_width = max(
        10.0,
        0.12
        * len(
            labels
        ),
    )

    figure, axis = plt.subplots(
        figsize=(
            figure_width,
            7.0,
        )
    )

    dendrogram(
        linkage_matrix,
        labels=labels,
        leaf_rotation=90,
        leaf_font_size=5,
        ax=axis,
        color_threshold=cutoff,
    )

    axis.axhline(
        cutoff,
        color="black",
        linestyle="--",
        linewidth=1.0,
    )

    axis.set_ylabel(
        distance_label
    )

    axis.set_title(
        title
    )

    figure.tight_layout()

    figure.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(
        figure
    )

# ============================================================
# CCP4 export
# ============================================================

def write_ccp4_map(
    density: np.ndarray,
    grid: dict,
    output_file: Path,
):
    """
    Write a regular Cartesian probability-density grid as an
    MRC/CCP4-compatible map.
    """
    if mrcfile is None:
        raise RuntimeError(
            "mrcfile is not available in this Python environment."
        )

    density = np.asarray(
        density,
        dtype=np.float32,
    )

    if density.ndim != 3:
        raise ValueError(
            "Density array must be three-dimensional."
        )

    if not np.all(
        np.isfinite(
            density
        )
    ):
        raise ValueError(
            "Density array contains non-finite values."
        )

    origin = np.asarray(
        grid["origin"],
        dtype=float,
    )

    spacing = float(
        grid["spacing"]
    )

    if origin.shape != (3,):
        raise ValueError(
            "Grid origin must contain exactly three coordinates."
        )

    if spacing <= 0.0:
        raise ValueError(
            "Grid spacing must be positive."
        )

    output_file = Path(
        output_file
    )

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    mrc_data = np.transpose(
        density,
        (
            2,
            1,
            0,
        ),
    )

    with mrcfile.new(
        str(
            output_file
        ),
        overwrite=True,
    ) as mrc:

        mrc.set_data(
            mrc_data.astype(
                np.float32,
                copy=False,
            )
        )

        mrc.voxel_size = (
            spacing,
            spacing,
            spacing,
        )

        mrc.header.origin.x = float(
            origin[0]
        )

        mrc.header.origin.y = float(
            origin[1]
        )

        mrc.header.origin.z = float(
            origin[2]
        )

        mrc.header.mapc = 1
        mrc.header.mapr = 2
        mrc.header.maps = 3

        mrc.update_header_stats()


def save_numpy_density(
    density: np.ndarray,
    output_file: Path,
):
    """
    Always write a lossless NumPy copy even if CCP4 export is unavailable.
    """
    np.save(
        output_file,
        density.astype(
            np.float32
        ),
    )


def save_grid_metadata(
    grid: dict,
    output_file: Path,
):
    data = {
        "origin": [
            float(value)
            for value in grid[
                "origin"
            ]
        ],
        "maximum": [
            float(value)
            for value in grid[
                "maximum"
            ]
        ],
        "spacing": float(
            grid[
                "spacing"
            ]
        ),
        "shape": [
            int(value)
            for value in grid[
                "shape"
            ]
        ],
    }

    output_file.write_text(
        json.dumps(
            data,
            indent=2,
        )
    )


# ============================================================
# Structure export
# ============================================================

def write_aligned_pose_archive(
    pose_records: list[dict],
    output_file: Path,
):
    with structure.StructureWriter(
        str(
            output_file
        )
    ) as writer:
        for record in pose_records:
            ligand_st = record[
                "structure"
            ]

            ligand_st.title = (
                f"{record['source_name']}_"
                f"pose{record['pose_number']:03d}_"
                f"family{record['global_family_id']:03d}"
            )

            ligand_st.property[
                "s_user_source_name"
            ] = record[
                "source_name"
            ]

            ligand_st.property[
                "i_user_pose_number"
            ] = int(
                record[
                    "pose_number"
                ]
            )

            ligand_st.property[
                "i_user_global_family"
            ] = int(
                record[
                    "global_family_id"
                ]
            )

            writer.append(
                ligand_st
            )

def residue_one_letter_code(
    residue_name: str,
) -> str:
    """
    Convert a prepared/PDB residue name to a one-letter code.

    Unknown protein residues are retained as X so that residue order is
    preserved, but X-X is not treated as an informative identity match.
    """
    residue_name = (
        residue_name
        or ""
    ).strip().upper()

    return THREE_TO_ONE.get(
        residue_name,
        "X",
    )


def extract_ca_chains(
    receptor_st,
    atom_name: str = "CA",
    minimum_chain_length: int = 25,
) -> dict[str, list[dict]]:
    """
    Extract ordered C-alpha residue records for every plausible protein
    chain.

    PDB chain identifiers are used only to group atoms within one
    structure. They are not required to agree between structures.
    """
    atom_name = atom_name.strip().upper()

    chains = defaultdict(
        list
    )

    seen_residues = set()

    for atom_order, atom in enumerate(
        receptor_st.atom
    ):
        current_atom_name = (
            getattr(
                atom,
                "pdbname",
                "",
            )
            or ""
        ).strip().upper()

        if current_atom_name != atom_name:
            continue

        chain_id = (
            getattr(
                atom,
                "chain",
                "",
            )
            or ""
        ).strip()

        residue_number = int(
            getattr(
                atom,
                "resnum",
                0,
            )
        )

        insertion_code = (
            getattr(
                atom,
                "inscode",
                "",
            )
            or ""
        ).strip()

        residue_name = (
            getattr(
                atom,
                "pdbres",
                "",
            )
            or ""
        ).strip().upper()

        residue_key = (
            chain_id,
            residue_number,
            insertion_code,
        )

        if residue_key in seen_residues:
            continue

        seen_residues.add(
            residue_key
        )

        chains[
            chain_id
        ].append(
            {
                "chain_id": chain_id,
                "resnum": residue_number,
                "inscode": insertion_code,
                "resname": residue_name,
                "one_letter": (
                    residue_one_letter_code(
                        residue_name
                    )
                ),
                "atom_index": int(
                    atom.index
                ),
                "atom_order": atom_order,
                "coordinate": np.asarray(
                    [
                        atom.x,
                        atom.y,
                        atom.z,
                    ],
                    dtype=float,
                ),
            }
        )

    # Atom iteration order is generally already polymer order. Sorting
    # by residue number is useful where normal numbering is available,
    # while atom_order keeps insertion/duplicate cases deterministic.
    output = {}

    for chain_id, residues in chains.items():
        residues = sorted(
            residues,
            key=lambda residue: (
                residue[
                    "resnum"
                ],
                residue[
                    "inscode"
                ],
                residue[
                    "atom_order"
                ],
            ),
        )

        if len(
            residues
        ) >= minimum_chain_length:
            output[
                chain_id
            ] = residues

    return output

def global_sequence_alignment(
    first_sequence: str,
    second_sequence: str,
    match_score: float = 2.0,
    mismatch_score: float = -1.0,
    gap_score: float = -2.0,
) -> dict:
    """
    Needleman-Wunsch global sequence alignment.

    Returns the aligned strings and index correspondences between
    ungapped sequence positions.
    """
    n_first = len(
        first_sequence
    )

    n_second = len(
        second_sequence
    )

    scores = np.zeros(
        (
            n_first + 1,
            n_second + 1,
        ),
        dtype=float,
    )

    traceback = np.zeros(
        (
            n_first + 1,
            n_second + 1,
        ),
        dtype=np.int8,
    )

    # Traceback codes:
    # 1 = diagonal
    # 2 = up
    # 3 = left

    for first_index in range(
        1,
        n_first + 1,
    ):
        scores[
            first_index,
            0,
        ] = (
            first_index
            * gap_score
        )

        traceback[
            first_index,
            0,
        ] = 2

    for second_index in range(
        1,
        n_second + 1,
    ):
        scores[
            0,
            second_index,
        ] = (
            second_index
            * gap_score
        )

        traceback[
            0,
            second_index,
        ] = 3

    for first_index in range(
        1,
        n_first + 1,
    ):
        first_letter = first_sequence[
            first_index - 1
        ]

        for second_index in range(
            1,
            n_second + 1,
        ):
            second_letter = second_sequence[
                second_index - 1
            ]

            is_informative_match = (
                first_letter == second_letter
                and first_letter != "X"
            )

            substitution_score = (
                match_score
                if is_informative_match
                else mismatch_score
            )

            diagonal = (
                scores[
                    first_index - 1,
                    second_index - 1,
                ]
                + substitution_score
            )

            up = (
                scores[
                    first_index - 1,
                    second_index,
                ]
                + gap_score
            )

            left = (
                scores[
                    first_index,
                    second_index - 1,
                ]
                + gap_score
            )

            candidates = (
                diagonal,
                up,
                left,
            )

            best_code = int(
                np.argmax(
                    candidates
                )
            ) + 1

            scores[
                first_index,
                second_index,
            ] = candidates[
                best_code - 1
            ]

            traceback[
                first_index,
                second_index,
            ] = best_code

    first_index = n_first
    second_index = n_second

    aligned_first = []
    aligned_second = []
    index_pairs = []

    while (
        first_index > 0
        or second_index > 0
    ):
        direction = traceback[
            first_index,
            second_index,
        ]

        if (
            first_index > 0
            and second_index > 0
            and direction == 1
        ):
            aligned_first.append(
                first_sequence[
                    first_index - 1
                ]
            )

            aligned_second.append(
                second_sequence[
                    second_index - 1
                ]
            )

            index_pairs.append(
                (
                    first_index - 1,
                    second_index - 1,
                )
            )

            first_index -= 1
            second_index -= 1

        elif (
            first_index > 0
            and (
                second_index == 0
                or direction == 2
            )
        ):
            aligned_first.append(
                first_sequence[
                    first_index - 1
                ]
            )

            aligned_second.append(
                "-"
            )

            first_index -= 1

        else:
            aligned_first.append(
                "-"
            )

            aligned_second.append(
                second_sequence[
                    second_index - 1
                ]
            )

            second_index -= 1

    aligned_first = "".join(
        reversed(
            aligned_first
        )
    )

    aligned_second = "".join(
        reversed(
            aligned_second
        )
    )

    index_pairs = list(
        reversed(
            index_pairs
        )
    )

    informative_pairs = [
        (
            first_position,
            second_position,
        )
        for (
            first_position,
            second_position,
        ) in index_pairs
        if (
            first_sequence[
                first_position
            ] != "X"
            and second_sequence[
                second_position
            ] != "X"
        )
    ]

    identical_pairs = [
        (
            first_position,
            second_position,
        )
        for (
            first_position,
            second_position,
        ) in informative_pairs
        if (
            first_sequence[
                first_position
            ]
            == second_sequence[
                second_position
            ]
        )
    ]

    aligned_pair_count = len(
        index_pairs
    )

    informative_pair_count = len(
        informative_pairs
    )

    identity = (
        len(
            identical_pairs
        )
        / informative_pair_count
        if informative_pair_count
        else 0.0
    )

    first_coverage = (
        aligned_pair_count
        / n_first
        if n_first
        else 0.0
    )

    second_coverage = (
        aligned_pair_count
        / n_second
        if n_second
        else 0.0
    )

    return {
        "score": float(
            scores[
                n_first,
                n_second,
            ]
        ),
        "aligned_first": aligned_first,
        "aligned_second": aligned_second,
        "index_pairs": index_pairs,
        "identity": float(
            identity
        ),
        "first_coverage": float(
            first_coverage
        ),
        "second_coverage": float(
            second_coverage
        ),
        "aligned_pair_count": int(
            aligned_pair_count
        ),
        "identical_pair_count": int(
            len(
                identical_pairs
            )
        ),
    }

def select_best_sequence_matched_chain_pair(
    reference_chains: dict[str, list[dict]],
    mobile_chains: dict[str, list[dict]],
    minimum_identity: float = 0.50,
    minimum_reference_coverage: float = 0.50,
) -> dict:
    """
    Compare every reference/mobile chain pair and choose the pair with
    the strongest sequence correspondence.

    Chain identifiers may differ completely between structures.
    """
    candidates = []

    for (
        reference_chain_id,
        reference_residues,
    ) in reference_chains.items():
        reference_sequence = "".join(
            residue[
                "one_letter"
            ]
            for residue in reference_residues
        )

        for (
            mobile_chain_id,
            mobile_residues,
        ) in mobile_chains.items():
            mobile_sequence = "".join(
                residue[
                    "one_letter"
                ]
                for residue in mobile_residues
            )

            alignment = global_sequence_alignment(
                reference_sequence,
                mobile_sequence,
            )

            # Prioritize sequence identity and coverage of the
            # reference receptor chain. This avoids selecting a short
            # matching fragment or a similarly sized unrelated chain.
            selection_score = (
                alignment[
                    "identity"
                ]
                * alignment[
                    "first_coverage"
                ]
                * math.sqrt(
                    max(
                        alignment[
                            "aligned_pair_count"
                        ],
                        1,
                    )
                )
            )

            candidates.append(
                {
                    "reference_chain_id": (
                        reference_chain_id
                    ),
                    "mobile_chain_id": (
                        mobile_chain_id
                    ),
                    "reference_residues": (
                        reference_residues
                    ),
                    "mobile_residues": (
                        mobile_residues
                    ),
                    "alignment": alignment,
                    "selection_score": float(
                        selection_score
                    ),
                }
            )

    if not candidates:
        raise RuntimeError(
            "No plausible protein chains containing C-alpha "
            "atoms were found."
        )

    candidates.sort(
        key=lambda candidate: (
            candidate[
                "selection_score"
            ],
            candidate[
                "alignment"
            ][
                "identity"
            ],
            candidate[
                "alignment"
            ][
                "aligned_pair_count"
            ],
        ),
        reverse=True,
    )

    best = candidates[
        0
    ]

    if (
        best[
            "alignment"
        ][
            "identity"
        ]
        < minimum_identity
    ):
        raise RuntimeError(
            "No receptor-chain pair passed the minimum "
            "sequence identity. Best chain pair "
            f"{best['reference_chain_id']!r} -> "
            f"{best['mobile_chain_id']!r} had identity "
            f"{best['alignment']['identity']:.3f}."
        )

    if (
        best[
            "alignment"
        ][
            "first_coverage"
        ]
        < minimum_reference_coverage
    ):
        raise RuntimeError(
            "No receptor-chain pair passed the minimum "
            "reference-chain sequence coverage. Best chain pair "
            f"{best['reference_chain_id']!r} -> "
            f"{best['mobile_chain_id']!r} had reference coverage "
            f"{best['alignment']['first_coverage']:.3f}."
        )

    # Record the runner-up margin as useful evidence of whether the
    # chain choice was unambiguous.
    if len(
        candidates
    ) > 1:
        runner_up_score = candidates[
            1
        ][
            "selection_score"
        ]

    else:
        runner_up_score = np.nan

    best[
        "runner_up_selection_score"
    ] = float(
        runner_up_score
    )

    return best

def sequence_mapped_pocket_coordinates(
    chain_match: dict,
    reference_ligand_coordinates: np.ndarray,
    pocket_cutoff: float,
) -> dict:
    """
    Select sequence-corresponding C-alpha pairs based on proximity of
    the reference-chain C-alpha to the reference ligand.

    The mobile ligand and mobile residue numbering are not used.
    """
    reference_residues = chain_match[
        "reference_residues"
    ]

    mobile_residues = chain_match[
        "mobile_residues"
    ]

    index_pairs = chain_match[
        "alignment"
    ][
        "index_pairs"
    ]

    reference_coordinates = []
    mobile_coordinates = []
    residue_pairs = []

    for (
        reference_index,
        mobile_index,
    ) in index_pairs:
        reference_residue = reference_residues[
            reference_index
        ]

        mobile_residue = mobile_residues[
            mobile_index
        ]

        reference_coordinate = (
            reference_residue[
                "coordinate"
            ]
        )

        minimum_ligand_distance = float(
            np.min(
                np.linalg.norm(
                    reference_ligand_coordinates
                    - reference_coordinate[
                        None,
                        :
                    ],
                    axis=1,
                )
            )
        )

        if (
            minimum_ligand_distance
            > pocket_cutoff
        ):
            continue

        reference_coordinates.append(
            reference_coordinate
        )

        mobile_coordinates.append(
            mobile_residue[
                "coordinate"
            ]
        )

        residue_pairs.append(
            {
                "reference_sequence_index": int(
                    reference_index
                ),
                "mobile_sequence_index": int(
                    mobile_index
                ),
                "reference_chain": (
                    reference_residue[
                        "chain_id"
                    ]
                ),
                "reference_resnum": (
                    reference_residue[
                        "resnum"
                    ]
                ),
                "reference_inscode": (
                    reference_residue[
                        "inscode"
                    ]
                ),
                "reference_resname": (
                    reference_residue[
                        "resname"
                    ]
                ),
                "mobile_chain": (
                    mobile_residue[
                        "chain_id"
                    ]
                ),
                "mobile_resnum": (
                    mobile_residue[
                        "resnum"
                    ]
                ),
                "mobile_inscode": (
                    mobile_residue[
                        "inscode"
                    ]
                ),
                "mobile_resname": (
                    mobile_residue[
                        "resname"
                    ]
                ),
                "reference_ligand_distance": (
                    minimum_ligand_distance
                ),
            }
        )

    return {
        "reference_coordinates": np.asarray(
            reference_coordinates,
            dtype=float,
        ),
        "mobile_coordinates": np.asarray(
            mobile_coordinates,
            dtype=float,
        ),
        "residue_pairs": residue_pairs,
    }

def robust_kabsch_transform(
    mobile_coordinates: np.ndarray,
    reference_coordinates: np.ndarray,
    minimum_atoms: int,
    outlier_cutoff: float = 3.0,
) -> dict:
    """
    Fit matched coordinates, remove gross residual outliers once, and
    refit.

    This does not repeatedly optimize the selected set and therefore
    remains easy to diagnose.
    """
    if (
        len(
            mobile_coordinates
        )
        < minimum_atoms
    ):
        raise RuntimeError(
            f"Only {len(mobile_coordinates)} matched pocket atoms "
            f"were available; at least {minimum_atoms} are required."
        )

    initial_rotation, initial_translation = (
        kabsch_transform(
            mobile_coordinates,
            reference_coordinates,
        )
    )

    initially_aligned = apply_transform(
        mobile_coordinates,
        initial_rotation,
        initial_translation,
    )

    initial_residuals = np.linalg.norm(
        initially_aligned
        - reference_coordinates,
        axis=1,
    )

    retained_mask = (
        initial_residuals
        <= outlier_cutoff
    )

    if int(
        retained_mask.sum()
    ) < minimum_atoms:
        raise RuntimeError(
            "Sequence correspondence was found, but only "
            f"{int(retained_mask.sum())} pocket atoms remained "
            f"within the {outlier_cutoff:.2f} A coordinate-outlier "
            f"cutoff; at least {minimum_atoms} are required."
        )

    retained_mobile = mobile_coordinates[
        retained_mask
    ]

    retained_reference = reference_coordinates[
        retained_mask
    ]

    final_rotation, final_translation = (
        kabsch_transform(
            retained_mobile,
            retained_reference,
        )
    )

    final_rmsd = alignment_rmsd(
        mobile_coordinates=retained_mobile,
        reference_coordinates=retained_reference,
        rotation=final_rotation,
        translation=final_translation,
    )

    final_aligned_all = apply_transform(
        mobile_coordinates,
        final_rotation,
        final_translation,
    )

    final_all_residuals = np.linalg.norm(
        final_aligned_all
        - reference_coordinates,
        axis=1,
    )

    return {
        "rotation": final_rotation,
        "translation": final_translation,
        "fit_rmsd": float(
            final_rmsd
        ),
        "retained_mask": retained_mask,
        "initial_residuals": initial_residuals,
        "final_all_residuals": final_all_residuals,
        "n_initial_atoms": int(
            len(
                mobile_coordinates
            )
        ),
        "n_retained_atoms": int(
            retained_mask.sum()
        ),
    }

def pairwise_centroid_distance_matrix(
    pose_records: list[dict],
) -> np.ndarray:
    """
    Calculate pairwise Euclidean distances between aligned ligand
    centroids.

    This provides the first-level campaign definition of spatial pose
    families and does not require atom correspondence between ligands.
    """
    if not pose_records:
        raise ValueError(
            "Cannot calculate centroid distances without poses."
        )

    centroids = np.asarray(
        [
            record[
                "centroid"
            ]
            for record in pose_records
        ],
        dtype=float,
    )

    if centroids.ndim != 2 or centroids.shape[1] != 3:
        raise ValueError(
            "Pose centroids must have shape (n_poses, 3)."
        )

    if len(
        centroids
    ) == 1:
        return np.zeros(
            (
                1,
                1,
            ),
            dtype=float,
        )

    return squareform(
        pdist(
            centroids,
            metric="euclidean",
        )
    )

def ligand_orientation_descriptors(
    coordinates: np.ndarray,
) -> dict:
    """
    Calculate translation-independent ligand shape and orientation
    descriptors from aligned heavy-atom coordinates.

    The first principal axis is unsigned: axis u and axis -u represent
    the same geometric direction. Contact fingerprints are therefore
    needed to distinguish many 180-degree flips.
    """
    coordinates = np.asarray(
        coordinates,
        dtype=float,
    )

    if (
        coordinates.ndim != 2
        or coordinates.shape[1] != 3
        or len(
            coordinates
        ) < 2
    ):
        raise ValueError(
            "Ligand coordinates must have shape (n, 3) "
            "with at least two atoms."
        )

    centroid = coordinates.mean(
        axis=0
    )

    centered = (
        coordinates
        - centroid
    )

    covariance = (
        centered.T
        @ centered
    ) / max(
        len(
            centered
        ) - 1,
        1,
    )

    eigenvalues, eigenvectors = np.linalg.eigh(
        covariance
    )

    order = np.argsort(
        eigenvalues
    )[
        ::-1
    ]

    eigenvalues = eigenvalues[
        order
    ]

    eigenvectors = eigenvectors[
        :,
        order
    ]

    principal_axis = eigenvectors[
        :,
        0
    ]

    principal_axis /= max(
        np.linalg.norm(
            principal_axis
        ),
        1.0e-12,
    )

    total_variance = float(
        eigenvalues.sum()
    )

    axis_anisotropy = (
        float(
            (
                eigenvalues[0]
                - eigenvalues[1]
            )
            / total_variance
        )
        if total_variance > 0.0
        else 0.0
    )

    radius_of_gyration = float(
        np.sqrt(
            np.mean(
                np.sum(
                    centered ** 2,
                    axis=1,
                )
            )
        )
    )

    return {
        "centroid": centroid,
        "principal_axis": principal_axis,
        "eigenvalues": eigenvalues,
        "axis_anisotropy": axis_anisotropy,
        "radius_of_gyration": radius_of_gyration,
    }


def unsigned_axis_angle_degrees(
    first_axis: np.ndarray,
    second_axis: np.ndarray,
) -> float:
    """
    Angle between unsigned principal axes in the range 0-90 degrees.
    """
    first_axis = np.asarray(
        first_axis,
        dtype=float,
    )

    second_axis = np.asarray(
        second_axis,
        dtype=float,
    )

    denominator = (
        np.linalg.norm(
            first_axis
        )
        * np.linalg.norm(
            second_axis
        )
    )

    if denominator <= 0.0:
        return 90.0

    cosine = float(
        np.dot(
            first_axis,
            second_axis,
        )
        / denominator
    )

    cosine = float(
        np.clip(
            abs(
                cosine
            ),
            0.0,
            1.0,
        )
    )

    return float(
        np.degrees(
            np.arccos(
                cosine
            )
        )
    )

def _normalise_element(
    atom,
) -> str:
    return (
        getattr(
            atom,
            "element",
            "",
        )
        or ""
    ).strip().upper()


def _residue_identity_from_atom(
    atom,
) -> tuple:
    return (
        (
            getattr(
                atom,
                "chain",
                "",
            )
            or ""
        ).strip(),
        int(
            getattr(
                atom,
                "resnum",
                0,
            )
        ),
        (
            getattr(
                atom,
                "inscode",
                "",
            )
            or ""
        ).strip(),
    )


def build_mobile_to_reference_position_map(
    chain_match: dict,
) -> dict[tuple, int]:
    """
    Map mobile receptor residue identities onto zero-based reference
    sequence positions using the selected sequence alignment.
    """
    mobile_residues = chain_match[
        "mobile_residues"
    ]

    mapping = {}

    for (
        reference_index,
        mobile_index,
    ) in chain_match[
        "alignment"
    ][
        "index_pairs"
    ]:
        mobile_residue = mobile_residues[
            mobile_index
        ]

        mobile_key = (
            mobile_residue[
                "chain_id"
            ],
            int(
                mobile_residue[
                    "resnum"
                ]
            ),
            mobile_residue[
                "inscode"
            ],
        )

        mapping[
            mobile_key
        ] = int(
            reference_index
        )

    return mapping


def receptor_contact_atoms(
    aligned_receptor,
    mobile_to_reference_position: dict,
) -> list[dict]:
    """
    Create a flat receptor-heavy-atom table carrying canonical
    reference-sequence positions.
    """
    records = []

    for atom in aligned_receptor.atom:
        element = _normalise_element(
            atom
        )

        if element == "H":
            continue

        mobile_key = _residue_identity_from_atom(
            atom
        )

        reference_position = (
            mobile_to_reference_position.get(
                mobile_key
            )
        )

        if reference_position is None:
            continue

        records.append(
            {
                "reference_position": int(
                    reference_position
                ),
                "element": element,
                "coordinate": np.asarray(
                    [
                        atom.x,
                        atom.y,
                        atom.z,
                    ],
                    dtype=float,
                ),
            }
        )

    return records


def pose_contact_fingerprint(
    ligand_coordinates: np.ndarray,
    ligand_elements: list[str],
    receptor_atoms: list[dict],
    contact_cutoff: float,
    polar_contact_cutoff: float,
) -> set[str]:
    """
    Return a sparse set of sequence-mapped receptor-contact features.

    Examples:
        R0123:contact
        R0123:polar_candidate
    """
    if not receptor_atoms:
        return set()

    receptor_coordinates = np.asarray(
        [
            record[
                "coordinate"
            ]
            for record in receptor_atoms
        ],
        dtype=float,
    )

    ligand_coordinates = np.asarray(
        ligand_coordinates,
        dtype=float,
    )

    distances = cdist(
        ligand_coordinates,
        receptor_coordinates,
    )

    ligand_elements = [
        str(
            element
        ).strip().upper()
        for element in ligand_elements
    ]

    features = set()

    ligand_indices, receptor_indices = np.where(
        distances <= contact_cutoff
    )

    polar_elements = {
        "N",
        "O",
        "S",
    }

    for ligand_index, receptor_index in zip(
        ligand_indices,
        receptor_indices,
    ):
        receptor_record = receptor_atoms[
            receptor_index
        ]

        reference_position = int(
            receptor_record[
                "reference_position"
            ]
        )

        feature_prefix = (
            f"R{reference_position + 1:04d}"
        )

        features.add(
            f"{feature_prefix}:contact"
        )

        distance = float(
            distances[
                ligand_index,
                receptor_index,
            ]
        )

        if (
            distance
            <= polar_contact_cutoff
            and ligand_elements[
                ligand_index
            ] in polar_elements
            and receptor_record[
                "element"
            ] in polar_elements
        ):
            features.add(
                f"{feature_prefix}:polar_candidate"
            )

    return features


def jaccard_distance(
    first_features: set[str],
    second_features: set[str],
) -> float:
    union = (
        first_features
        | second_features
    )

    if not union:
        return 0.0

    intersection = (
        first_features
        & second_features
    )

    return float(
        1.0
        - len(
            intersection
        )
        / len(
            union
        )
    )

def second_level_distance_matrix(
    pose_records: list[dict],
    orientation_weight: float,
    contact_weight: float,
    minimum_axis_anisotropy: float,
) -> np.ndarray:
    """
    Combine orientation and contact dissimilarity within one spatial
    family.

    If either ligand is too isotropic for a stable principal axis,
    contact distance receives the full weight.
    """
    n_poses = len(
        pose_records
    )

    matrix = np.zeros(
        (
            n_poses,
            n_poses,
        ),
        dtype=float,
    )

    total_weight = (
        orientation_weight
        + contact_weight
    )

    if total_weight <= 0.0:
        raise ValueError(
            "Second-level clustering weights must sum "
            "to a positive value."
        )

    orientation_weight /= total_weight
    contact_weight /= total_weight

    for first_index in range(
        n_poses
    ):
        first = pose_records[
            first_index
        ]

        for second_index in range(
            first_index + 1,
            n_poses,
        ):
            second = pose_records[
                second_index
            ]

            contact_distance = jaccard_distance(
                first[
                    "contact_fingerprint"
                ],
                second[
                    "contact_fingerprint"
                ],
            )

            axes_are_stable = (
                first[
                    "axis_anisotropy"
                ]
                >= minimum_axis_anisotropy
                and second[
                    "axis_anisotropy"
                ]
                >= minimum_axis_anisotropy
            )

            if axes_are_stable:
                angle = unsigned_axis_angle_degrees(
                    first[
                        "principal_axis"
                    ],
                    second[
                        "principal_axis"
                    ],
                )

                orientation_distance = (
                    angle
                    / 90.0
                )

                combined_distance = (
                    orientation_weight
                    * orientation_distance
                    + contact_weight
                    * contact_distance
                )

            else:
                combined_distance = (
                    contact_distance
                )

            matrix[
                first_index,
                second_index,
            ] = combined_distance

            matrix[
                second_index,
                first_index,
            ] = combined_distance

    return matrix


def assign_second_level_subfamilies(
    all_pose_records: list[dict],
    distance_cutoff: float,
    orientation_weight: float,
    contact_weight: float,
    minimum_axis_anisotropy: float,
    linkage_method: str = "average",
) -> dict[int, np.ndarray]:
    """
    Cluster poses separately inside each global spatial family.

    Subfamily IDs are local to their parent spatial family.
    """
    matrices = {}

    spatial_family_ids = sorted(
        {
            int(
                record[
                    "global_family_id"
                ]
            )
            for record in all_pose_records
        }
    )

    for spatial_family_id in spatial_family_ids:
        family_indices = [
            index
            for index, record in enumerate(
                all_pose_records
            )
            if int(
                record[
                    "global_family_id"
                ]
            )
            == spatial_family_id
        ]

        family_records = [
            all_pose_records[
                index
            ]
            for index in family_indices
        ]

        if len(
            family_records
        ) == 1:
            all_pose_records[
                family_indices[0]
            ][
                "orientation_contact_subfamily_id"
            ] = 1

            matrices[
                spatial_family_id
            ] = np.zeros(
                (
                    1,
                    1,
                ),
                dtype=float,
            )

            continue

        distance_matrix = (
            second_level_distance_matrix(
                pose_records=family_records,
                orientation_weight=orientation_weight,
                contact_weight=contact_weight,
                minimum_axis_anisotropy=(
                    minimum_axis_anisotropy
                ),
            )
        )

        linkage_matrix = linkage(
            squareform(
                distance_matrix,
                checks=False,
            ),
            method=linkage_method,
        )

        subfamily_ids = fcluster(
            linkage_matrix,
            t=distance_cutoff,
            criterion="distance",
        )

        # Renumber by population.
        unique_ids, counts = np.unique(
            subfamily_ids,
            return_counts=True,
        )

        order = unique_ids[
            np.argsort(
                -counts
            )
        ]

        remapping = {
            int(
                old_id
            ): new_id
            for new_id, old_id in enumerate(
                order,
                start=1,
            )
        }

        for global_index, raw_id in zip(
            family_indices,
            subfamily_ids,
        ):
            all_pose_records[
                global_index
            ][
                "orientation_contact_subfamily_id"
            ] = remapping[
                int(
                    raw_id
                )
            ]

        matrices[
            spatial_family_id
        ] = distance_matrix

    return matrices


def export_second_level_pymol(
    pose_records: list[dict],
    aligned_receptor_file: Path,
    output_dir: Path,
    alignment_reference_name: str | None = None,
) -> Path:
    """
    Export the second-level GlideEM landscape for PyMOL.

    Output:
        - one reference-receptor PDB;
        - one ligand PDB per GlideEM pose;
        - one PyMOL script grouping poses as:

              spatial family
                  orientation/contact subfamily
                      individual source pose

    One structure is written per PDB because Schrodinger's PDB writer
    does not support writing multiple Structure objects to one file.
    """
    output_dir = (
        Path(
            output_dir
        )
        .expanduser()
        .resolve()
    )

    aligned_receptor_file = (
        Path(
            aligned_receptor_file
        )
        .expanduser()
        .resolve()
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not aligned_receptor_file.is_file():
        raise FileNotFoundError(
            "Aligned receptor archive was not found: "
            f"{aligned_receptor_file}"
        )

    if not pose_records:
        raise ValueError(
            "No pose records were supplied for "
            "second-level PyMOL export."
        )

    # --------------------------------------------------
    # Remove outputs from earlier exporter versions
    # --------------------------------------------------

    for pattern in (
        "*.pdb",
        "*.maegz",
        "*.pml",
    ):
        for stale_file in output_dir.glob(
            pattern
        ):
            stale_file.unlink()

    # --------------------------------------------------
    # Validate and group pose records
    # --------------------------------------------------

    required_pose_fields = {
        "source_name",
        "pose_number",
        "global_family_id",
        "orientation_contact_subfamily_id",
        "structure",
    }

    grouped_records = defaultdict(
        list
    )

    for record_index, record in enumerate(
        pose_records
    ):
        missing_fields = (
            required_pose_fields
            - set(
                record
            )
        )

        if missing_fields:
            raise RuntimeError(
                f"Pose record {record_index} is missing fields: "
                f"{sorted(missing_fields)}"
            )

        spatial_family_id = int(
            record[
                "global_family_id"
            ]
        )

        subfamily_id = int(
            record[
                "orientation_contact_subfamily_id"
            ]
        )

        grouped_records[
            (
                spatial_family_id,
                subfamily_id,
            )
        ].append(
            record
        )

    if not grouped_records:
        raise RuntimeError(
            "No second-level pose groups were generated."
        )

    # --------------------------------------------------
    # Select one aligned receptor as structural context
    # --------------------------------------------------

    reference_receptor = None
    available_receptor_names = []

    for receptor_st in structure.StructureReader(
        str(
            aligned_receptor_file
        )
    ):
        source_name = str(
            receptor_st.property.get(
                "s_user_source_name",
                "",
            )
        )

        available_receptor_names.append(
            source_name
        )

        if (
            alignment_reference_name is not None
            and source_name
            == alignment_reference_name
        ):
            reference_receptor = receptor_st
            break

        if reference_receptor is None:
            # Fallback used only if no explicit reference name was
            # supplied or its metadata cannot be found.
            reference_receptor = receptor_st

    if reference_receptor is None:
        raise RuntimeError(
            "No receptor structures were found in "
            f"{aligned_receptor_file}"
        )

    if (
        alignment_reference_name is not None
        and str(
            reference_receptor.property.get(
                "s_user_source_name",
                "",
            )
        )
        != alignment_reference_name
    ):
        raise RuntimeError(
            "Could not find the requested alignment-reference "
            f"receptor '{alignment_reference_name}' in "
            f"{aligned_receptor_file}. Available source names: "
            f"{available_receptor_names}"
        )

    receptor_source_name = str(
        reference_receptor.property.get(
            "s_user_source_name",
            "reference",
        )
    )

    receptor_object_name = (
        "aligned_reference_receptor"
    )

    receptor_pdb = (
        output_dir
        / (
            f"{_safe_object_name(receptor_source_name)}"
            "_aligned_reference_receptor.pdb"
        )
    )

    # One Structure object, one PDB file.
    with structure.StructureWriter(
        str(
            receptor_pdb
        )
    ) as writer:
        writer.append(
            reference_receptor
        )

    if (
        not receptor_pdb.is_file()
        or receptor_pdb.stat().st_size == 0
    ):
        raise RuntimeError(
            "Failed to export the aligned reference receptor: "
            f"{receptor_pdb}"
        )

    # --------------------------------------------------
    # Initialise the PyMOL script
    # --------------------------------------------------

    script_lines = [
        "reinitialize",
        "bg_color white",
        "",
        (
            f"load {receptor_pdb.name}, "
            "receptor"
        ),
        "hide everything, receptor",
        "show cartoon, receptor",
        "color gray80, receptor",
        "",
    ]

    # --------------------------------------------------
    # Export one PDB per ligand pose
    # --------------------------------------------------

    pose_object_names = []
    exported_pose_count = 0

    for group_index, key in enumerate(
        sorted(
            grouped_records
        )
    ):
        (
            spatial_family_id,
            subfamily_id,
        ) = key

        spatial_group_name = (
            f"spatial_family_"
            f"{spatial_family_id:03d}"
        )

        subfamily_group_name = (
            f"sf{spatial_family_id:03d}"
            f"_sub{subfamily_id:03d}"
        )

        colour = PYMOL_COLOURS[
            group_index
            % len(
                PYMOL_COLOURS
            )
        ]

        subfamily_objects = []

        group_records = sorted(
            grouped_records[
                key
            ],
            key=lambda record: (
                str(
                    record[
                        "source_name"
                    ]
                ),
                int(
                    record[
                        "pose_number"
                    ]
                ),
            ),
        )

        for record in group_records:
            source_name = _safe_object_name(
                str(
                    record[
                        "source_name"
                    ]
                )
            )

            pose_number = int(
                record[
                    "pose_number"
                ]
            )

            object_name = (
                f"pose__"
                f"{source_name}_"
                f"{pose_number:03d}_"
                f"sf{spatial_family_id:03d}_"
                f"sub{subfamily_id:03d}"
            )

            pdb_file = (
                output_dir
                / f"{object_name}.pdb"
            )

            ligand_st = record[
                "structure"
            ].copy()

            ligand_st.title = (
                object_name
            )

            with structure.StructureWriter(
                str(
                    pdb_file
                )
            ) as writer:
                writer.append(
                    ligand_st
                )

            if (
                not pdb_file.is_file()
                or pdb_file.stat().st_size == 0
            ):
                raise RuntimeError(
                    "Failed to export ligand pose PDB: "
                    f"{pdb_file}"
                )

            script_lines.extend(
                [
                    "",
                    (
                        f"load {pdb_file.name}, "
                        f"{object_name}"
                    ),
                    (
                        f"show sticks, "
                        f"{object_name}"
                    ),
                    (
                        f"color {colour}, "
                        f"({object_name} and elem C)"
                    ),
                    (
                        f"color blue, "
                        f"({object_name} and elem N)"
                    ),
                    (
                        f"color red, "
                        f"({object_name} and elem O)"
                    ),
                    (
                        f"color yellow, "
                        f"({object_name} and elem S)"
                    ),
                    (
                        f"color orange, "
                        f"({object_name} and elem P)"
                    ),
                    (
                        f"color white, "
                        f"({object_name} and elem H)"
                    ),
                    (
                        f"color palecyan, "
                        f"({object_name} and elem F)"
                    ),
                    (
                        f"color green, "
                        f"({object_name} and elem Cl)"
                    ),
                    (
                        f"color brown, "
                        f"({object_name} and elem Br)"
                    ),
                    (
                        f"color purple, "
                        f"({object_name} and elem I)"
                    ),
                ]
            )

            subfamily_objects.append(
                object_name
            )

            exported_pose_count += 1

        # Group all individual pose objects belonging to this
        # orientation/contact subfamily in one operation.
        script_lines.append(
            (
                f"group {subfamily_group_name}, "
                + " ".join(
                    subfamily_objects
                )
            )
        )

        script_lines.append(
            (
                f"group {spatial_group_name}, "
                f"{subfamily_group_name}"
            )
        )

    if exported_pose_count != len(
        pose_records
    ):
        raise RuntimeError(
            "Second-level PyMOL export wrote an unexpected "
            f"number of poses: wrote {exported_pose_count}, "
            f"expected {len(pose_records)}."
        )

    # --------------------------------------------------
    # Whole-ensemble selection
    # --------------------------------------------------

    if pose_object_names:
        glideem_selection = (
            " or ".join(
                pose_object_names
            )
        )

        script_lines.extend(
            [
                "",
                "# Selection containing every aligned GlideEM pose.",
                (
                    "select glideem_pose_ensemble, "
                    f"{glideem_selection}"
                ),
                "",
                "# Protein residues surrounding the pose ensemble.",
                (
                    "select binding_site, "
                    "byres (aligned_reference_receptor "
                    "within 5.0 of glideem_pose_ensemble)"
                ),
                "show sticks, binding_site",
                "color gray80, (binding_site and elem C)",
                "color blue, (binding_site and elem N)",
                "color red, (binding_site and elem O)",
                "color yellow, (binding_site and elem S)",
                "color orange, (binding_site and elem P)",
                "color white, (binding_site and elem H)",
                "",
                "# Consistent thin-stick representation.",
                "set stick_radius, 0.15, glideem_pose_ensemble",
                "set stick_radius, 0.15, binding_site",
                "",
                "# Initial view centred on the campaign pose landscape.",
                "orient glideem_pose_ensemble",
                "zoom glideem_pose_ensemble, 8",
            ]
        )

    script_file = (
        output_dir
        / "inspect_second_level_landscape.pml"
    )

    script_file.write_text(
        "\n".join(
            script_lines
        )
        + "\n"
    )

    return script_file

def inferred_alternative_group(
    source_name: str,
) -> str:
    """
    Group explicit alternative labels such as:
        structure_52555_1
        structure_52555_2
        structure_51292
        structure_51292_alt1

    This is intentionally conservative and only handles the naming
    patterns currently used by the workflow.
    """
    if source_name.endswith(
        "_alt1"
    ):
        return source_name[
            :-5
        ]

    match = re.match(
        r"^(.*)_([12])$",
        source_name,
    )

    if match:
        return match.group(
            1
        )

    return source_name

def save_second_level_family_heatmap(
    pose_table: pd.DataFrame,
    output_file: Path,
):
    """
    Plot the fraction of each structure's complete GlideEM ensemble
    assigned to each orientation/contact subfamily.

    Columns are hierarchical:
        F1.1 = spatial family 1, subfamily 1
        F1.2 = spatial family 1, subfamily 2
        F2.1 = spatial family 2, subfamily 1
        ...

    Fractions are calculated relative to all poses from the structure,
    not merely poses within the parent spatial family.
    """
    data = pose_table.copy()

    data[
        "family_subfamily"
    ] = [
        (
            f"F{int(family_id)}."
            f"{int(subfamily_id)}"
        )
        for family_id, subfamily_id in zip(
            data[
                "global_spatial_family_id"
            ],
            data[
                "orientation_contact_subfamily_id"
            ],
        )
    ]

    counts = (
        data
        .groupby(
            [
                "source_name",
                "family_subfamily",
            ]
        )
        .size()
        .rename(
            "n_poses"
        )
        .reset_index()
    )

    structure_totals = (
        data
        .groupby(
            "source_name"
        )
        .size()
    )

    counts[
        "pose_fraction"
    ] = [
        (
            row.n_poses
            / structure_totals.loc[
                row.source_name
            ]
        )
        for row in counts.itertuples()
    ]

    matrix = (
        counts
        .pivot(
            index="source_name",
            columns="family_subfamily",
            values="pose_fraction",
        )
        .fillna(
            0.0
        )
    )

    # Sort hierarchically rather than lexically.
    def family_sort_key(
        label: str,
    ):
        family_part, subfamily_part = (
            label[
                1:
            ].split(
                "."
            )
        )

        return (
            int(
                family_part
            ),
            int(
                subfamily_part
            ),
        )

    ordered_columns = sorted(
        matrix.columns,
        key=family_sort_key,
    )

    matrix = matrix[
        ordered_columns
    ]

    # Keep the same structure ordering as the first-level plot:
    # most heterogeneous structures first.
    effective_counts = (
        data
        .groupby(
            "source_name"
        )[
            "family_subfamily"
        ]
        .value_counts(
            normalize=True
        )
        .groupby(
            level=0
        )
        .apply(
            lambda probabilities: float(
                np.exp(
                    -np.sum(
                        probabilities.to_numpy(
                            dtype=float
                        )
                        * np.log(
                            probabilities.to_numpy(
                                dtype=float
                            )
                        )
                    )
                )
            )
        )
    )

    structure_order = (
        effective_counts
        .sort_values(
            ascending=False
        )
        .index
        .tolist()
    )

    matrix = matrix.loc[
        structure_order
    ]

    matrix_values = matrix.to_numpy(
        dtype=float
    )

    colour_maximum = max(
        float(
            np.nanmax(
                matrix_values
            )
        ),
        0.10,
    )

    figure_width = min(
        28.0,
        max(
            10.0,
            0.42
            * matrix.shape[1]
            + 5.0,
        ),
    )

    figure_height = max(
        5.0,
        0.55
        * matrix.shape[0]
        + 2.0,
    )

    figure, axis = plt.subplots(
        figsize=(
            figure_width,
            figure_height,
        )
    )

    image = axis.imshow(
        matrix_values,
        aspect="auto",
        cmap="Blues",
        vmin=0.0,
        vmax=colour_maximum,
    )

    axis.set_xticks(
        np.arange(
            matrix.shape[1]
        )
    )

    axis.set_xticklabels(
        matrix.columns,
        rotation=90,
    )

    y_labels = [
        (
            f"{source_name} "
            f"(N_eff={effective_counts.loc[source_name]:.2f})"
        )
        for source_name in matrix.index
    ]

    axis.set_yticks(
        np.arange(
            matrix.shape[0]
        )
    )

    axis.set_yticklabels(
        y_labels
    )

    axis.set_xlabel(
        "Spatial family . orientation/contact subfamily"
    )

    axis.set_ylabel(
        "Cryo-EM structure / pose interpretation"
    )

    axis.set_title(
        "GlideEM second-level binding-mode landscape"
    )

    colour_bar = figure.colorbar(
        image,
        ax=axis,
    )

    colour_bar.set_label(
        "Fraction of structure's GlideEM poses"
    )

    # Draw boundaries between spatial families.
    previous_family = None

    for column_index, label in enumerate(
        matrix.columns
    ):
        family = int(
            label[
                1:
            ].split(
                "."
            )[
                0
            ]
        )

        if (
            previous_family is not None
            and family != previous_family
        ):
            axis.axvline(
                column_index - 0.5,
                color="black",
                linewidth=1.2,
            )

        previous_family = family

    figure.tight_layout()

    figure.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(
        figure
    )

def second_level_contact_summary(
    pose_table: pd.DataFrame,
) -> pd.DataFrame:
    """
    Summarize receptor contact-feature prevalence within every
    spatial/orientation-contact subfamily.
    """
    rows = []

    grouped = pose_table.groupby(
        [
            "global_spatial_family_id",
            "orientation_contact_subfamily_id",
        ]
    )

    for (
        spatial_family_id,
        subfamily_id,
    ), group in grouped:

        n_poses = len(
            group
        )

        feature_counts = defaultdict(
            int
        )

        for fingerprint_text in group[
            "contact_fingerprint"
        ].fillna(
            ""
        ):
            features = {
                feature
                for feature in str(
                    fingerprint_text
                ).split(
                    ";"
                )
                if feature
            }

            for feature in features:
                feature_counts[
                    feature
                ] += 1

        for feature, count in sorted(
            feature_counts.items(),
            key=lambda item: (
                -item[1],
                item[0],
            ),
        ):
            rows.append(
                {
                    "global_spatial_family_id": int(
                        spatial_family_id
                    ),
                    "orientation_contact_subfamily_id": int(
                        subfamily_id
                    ),
                    "subfamily_label": (
                        f"F{int(spatial_family_id)}."
                        f"{int(subfamily_id)}"
                    ),
                    "n_subfamily_poses": int(
                        n_poses
                    ),
                    "contact_feature": feature,
                    "n_poses_with_feature": int(
                        count
                    ),
                    "feature_fraction": float(
                        count
                        / n_poses
                    ),
                }
            )

    return pd.DataFrame(
        rows
    )