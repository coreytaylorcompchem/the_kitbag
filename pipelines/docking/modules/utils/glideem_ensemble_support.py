from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from scipy.spatial.distance import jensenshannon

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

# ============================================================
# Input preparation
# ============================================================

def _validate_pose_table(
    pose_table: pd.DataFrame,
):
    required_columns = {
        "source_name",
        "pose_number",
        "local_cluster_id",
        "global_spatial_family_id",
        "orientation_contact_subfamily_id",
    }

    missing_columns = (
        required_columns
        - set(
            pose_table.columns
        )
    )

    if missing_columns:
        raise ValueError(
            "Ensemble-support analysis requires columns: "
            f"{sorted(required_columns)}. Missing: "
            f"{sorted(missing_columns)}"
        )

    if pose_table.empty:
        raise ValueError(
            "Pose table is empty."
        )


def _spatial_label(
    family_id: int,
) -> str:
    return f"F{int(family_id)}"


def _joint_label(
    family_id: int,
    subfamily_id: int,
) -> str:
    return (
        f"F{int(family_id)}."
        f"{int(subfamily_id)}"
    )


def _mode_labels(
    pose_table: pd.DataFrame,
    mode_level: str,
) -> list[str]:
    if mode_level == "spatial":
        labels = sorted(
            {
                _spatial_label(
                    family_id
                )
                for family_id in pose_table[
                    "global_spatial_family_id"
                ].astype(
                    int
                )
            },
            key=lambda label: int(
                label[
                    1:
                ]
            ),
        )

    elif mode_level == "joint":
        pairs = sorted(
            {
                (
                    int(
                        family_id
                    ),
                    int(
                        subfamily_id
                    ),
                )
                for family_id, subfamily_id in zip(
                    pose_table[
                        "global_spatial_family_id"
                    ],
                    pose_table[
                        "orientation_contact_subfamily_id"
                    ],
                )
            }
        )

        labels = [
            _joint_label(
                family_id,
                subfamily_id,
            )
            for (
                family_id,
                subfamily_id,
            ) in pairs
        ]

    else:
        raise ValueError(
            f"Unknown mode level: {mode_level}"
        )

    return labels


def _pose_mode_labels(
    pose_table: pd.DataFrame,
    mode_level: str,
) -> pd.Series:
    if mode_level == "spatial":
        return pose_table[
            "global_spatial_family_id"
        ].astype(
            int
        ).map(
            _spatial_label
        )

    if mode_level == "joint":
        return pd.Series(
            [
                _joint_label(
                    family_id,
                    subfamily_id,
                )
                for family_id, subfamily_id in zip(
                    pose_table[
                        "global_spatial_family_id"
                    ],
                    pose_table[
                        "orientation_contact_subfamily_id"
                    ],
                )
            ],
            index=pose_table.index,
        )

    raise ValueError(
        f"Unknown mode level: {mode_level}"
    )


# ============================================================
# Effective sample sizes
# ============================================================

def inverse_simpson_effective_size(
    local_cluster_ids: pd.Series,
) -> float:
    """
    Effective number of local pose clusters:

        N_eff = 1 / sum(q_c ** 2)

    where q_c is the fraction of poses assigned to local cluster c.
    """
    valid_ids = (
        local_cluster_ids
        .dropna()
        .astype(
            str
        )
    )

    if valid_ids.empty:
        return 1.0

    probabilities = (
        valid_ids
        .value_counts(
            normalize=True
        )
        .to_numpy(
            dtype=float
        )
    )

    denominator = float(
        np.sum(
            probabilities ** 2
        )
    )

    if denominator <= 0.0:
        return 1.0

    return float(
        1.0
        / denominator
    )


def structure_effective_sizes(
    pose_table: pd.DataFrame,
    enabled: bool,
    minimum: float,
) -> pd.DataFrame:
    rows = []

    for source_name, group in pose_table.groupby(
        "source_name",
        sort=True,
    ):
        raw_nposes = int(
            len(
                group
            )
        )

        if enabled:
            effective_n = (
                inverse_simpson_effective_size(
                    group[
                        "local_cluster_id"
                    ]
                )
            )

        else:
            effective_n = float(
                raw_nposes
            )

        effective_n = max(
            float(
                minimum
            ),
            min(
                float(
                    raw_nposes
                ),
                float(
                    effective_n
                ),
            ),
        )

        rows.append(
            {
                "source_name": source_name,
                "raw_nposes": raw_nposes,
                "effective_nposes": effective_n,
                "redundancy_ratio": float(
                    effective_n
                    / raw_nposes
                    if raw_nposes
                    else np.nan
                ),
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================
# Dirichlet posteriors
# ============================================================

def observed_mode_proportions(
    pose_table: pd.DataFrame,
    mode_level: str,
    source_names: list[str],
    mode_labels: list[str],
) -> np.ndarray:
    working = pose_table.copy()

    working[
        "_mode_label"
    ] = _pose_mode_labels(
        working,
        mode_level,
    )

    matrix = np.zeros(
        (
            len(
                source_names
            ),
            len(
                mode_labels
            ),
        ),
        dtype=float,
    )

    mode_to_index = {
        label: index
        for index, label in enumerate(
            mode_labels
        )
    }

    source_to_index = {
        source_name: index
        for index, source_name in enumerate(
            source_names
        )
    }

    counts = (
        working
        .groupby(
            [
                "source_name",
                "_mode_label",
            ]
        )
        .size()
    )

    source_totals = (
        working
        .groupby(
            "source_name"
        )
        .size()
    )

    for (
        source_name,
        mode_label,
    ), count in counts.items():
        matrix[
            source_to_index[
                source_name
            ],
            mode_to_index[
                mode_label
            ],
        ] = (
            float(
                count
            )
            / float(
                source_totals.loc[
                    source_name
                ]
            )
        )

    return matrix


def draw_structure_posteriors(
    observed_proportions: np.ndarray,
    effective_sizes: np.ndarray,
    n_draws: int,
    alpha: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Output shape:

        n_draws x n_structures x n_modes
    """
    n_structures, n_modes = (
        observed_proportions.shape
    )

    draws = np.zeros(
        (
            n_draws,
            n_structures,
            n_modes,
        ),
        dtype=np.float32,
    )

    for structure_index in range(
        n_structures
    ):
        effective_counts = (
            observed_proportions[
                structure_index
            ]
            * effective_sizes[
                structure_index
            ]
        )

        parameters = (
            effective_counts
            + alpha
        )

        draws[
            :,
            structure_index,
            :,
        ] = rng.dirichlet(
            parameters,
            size=n_draws,
        ).astype(
            np.float32
        )

    return draws


# ============================================================
# Data-driven landscape similarity and weights
# ============================================================

def js_distance_matrix(
    probabilities: np.ndarray,
) -> np.ndarray:
    """
    probabilities shape:

        n_structures x n_modes
    """
    n_structures = probabilities.shape[
        0
    ]

    matrix = np.zeros(
        (
            n_structures,
            n_structures,
        ),
        dtype=float,
    )

    for first_index in range(
        n_structures
    ):
        for second_index in range(
            first_index + 1,
            n_structures,
        ):
            distance = float(
                jensenshannon(
                    probabilities[
                        first_index
                    ],
                    probabilities[
                        second_index
                    ],
                    base=2.0,
                )
            )

            if not np.isfinite(
                distance
            ):
                distance = 0.0

            matrix[
                first_index,
                second_index,
            ] = distance

            matrix[
                second_index,
                first_index,
            ] = distance

    return matrix


def positive_distance_median(
    distance_matrix: np.ndarray,
) -> float:
    values = distance_matrix[
        np.triu_indices_from(
            distance_matrix,
            k=1,
        )
    ]

    values = values[
        np.isfinite(
            values
        )
        & (
            values > 0.0
        )
    ]

    if values.size == 0:
        return 1.0

    return float(
        np.median(
            values
        )
    )


def similarity_kernel_weights(
    distance_matrix: np.ndarray,
    bandwidth_multiplier: float,
) -> tuple[np.ndarray, float]:
    """
    Continuous redundancy adjustment:

        S_st = exp(-D_st**2 / (2 h**2))
        raw_weight_s = 1 / sum_t S_st
    """
    base_bandwidth = positive_distance_median(
        distance_matrix
    )

    bandwidth = max(
        float(
            base_bandwidth
            * bandwidth_multiplier
        ),
        1.0e-8,
    )

    similarities = np.exp(
        -(
            distance_matrix ** 2
        )
        / (
            2.0
            * bandwidth ** 2
        )
    )

    raw_weights = (
        1.0
        / np.sum(
            similarities,
            axis=1,
        )
    )

    total = float(
        raw_weights.sum()
    )

    if total <= 0.0:
        weights = np.full(
            len(
                raw_weights
            ),
            1.0
            / len(
                raw_weights
            ),
            dtype=float,
        )

    else:
        weights = (
            raw_weights
            / total
        )

    return (
        weights,
        bandwidth,
    )


def campaign_support_draws(
    structure_draws: np.ndarray,
    bandwidth_multipliers: list[float],
) -> dict:
    """
    Similarity is learned from the complete joint-mode landscape for
    each draw.

    Returns equal-structure support and one adjusted result for each
    bandwidth multiplier.
    """
    n_draws, n_structures, n_modes = (
        structure_draws.shape
    )

    equal_support = np.mean(
        structure_draws,
        axis=1,
    ).astype(
        np.float32
    )

    result = {
        "equal": equal_support,
        "adjusted": {},
        "weights": {},
        "bandwidths": {},
        "mean_js_matrix": np.zeros(
            (
                n_structures,
                n_structures,
            ),
            dtype=float,
        ),
    }

    for multiplier in bandwidth_multipliers:
        result[
            "adjusted"
        ][
            multiplier
        ] = np.zeros(
            (
                n_draws,
                n_modes,
            ),
            dtype=np.float32,
        )

        result[
            "weights"
        ][
            multiplier
        ] = np.zeros(
            (
                n_draws,
                n_structures,
            ),
            dtype=np.float32,
        )

        result[
            "bandwidths"
        ][
            multiplier
        ] = np.zeros(
            n_draws,
            dtype=np.float32,
        )

    for draw_index in range(
        n_draws
    ):
        probabilities = structure_draws[
            draw_index
        ]

        distance_matrix = js_distance_matrix(
            probabilities
        )

        result[
            "mean_js_matrix"
        ] += distance_matrix

        for multiplier in bandwidth_multipliers:
            weights, bandwidth = (
                similarity_kernel_weights(
                    distance_matrix,
                    bandwidth_multiplier=multiplier,
                )
            )

            result[
                "weights"
            ][
                multiplier
            ][
                draw_index
            ] = weights

            result[
                "bandwidths"
            ][
                multiplier
            ][
                draw_index
            ] = bandwidth

            result[
                "adjusted"
            ][
                multiplier
            ][
                draw_index
            ] = (
                weights
                @ probabilities
            )

    result[
        "mean_js_matrix"
    ] /= float(
        n_draws
    )

    return result


def apply_learned_weights(
    structure_draws: np.ndarray,
    weights_by_multiplier: dict,
) -> dict:
    """
    Apply weights learned from the joint landscape to another mode
    level, such as spatial-family support.
    """
    results = {}

    for multiplier, weights in (
        weights_by_multiplier.items()
    ):
        results[
            multiplier
        ] = np.einsum(
            "bs,bsm->bm",
            weights,
            structure_draws,
        ).astype(
            np.float32
        )

    return results


# ============================================================
# Posterior summaries
# ============================================================

def posterior_summary(
    draws: np.ndarray,
    mode_labels: list[str],
    observed_proportions: np.ndarray,
    support_threshold: float,
    analysis_name: str,
) -> pd.DataFrame:
    """
    draws shape:

        n_draws x n_modes
    """
    n_draws, n_modes = draws.shape

    order_per_draw = np.argsort(
        -draws,
        axis=1,
    )

    rank_matrix = np.empty_like(
        order_per_draw
    )

    row_indices = np.arange(
        n_draws
    )[
        :,
        None
    ]

    rank_matrix[
        row_indices,
        order_per_draw
    ] = (
        np.arange(
            n_modes
        )[
            None,
            :
        ]
        + 1
    )

    rows = []

    for mode_index, mode_label in enumerate(
        mode_labels
    ):
        values = draws[
            :,
            mode_index
        ]

        observed_supporting_structures = int(
            np.sum(
                observed_proportions[
                    :,
                    mode_index
                ]
                >= support_threshold
            )
        )

        rows.append(
            {
                "analysis": analysis_name,
                "mode": mode_label,
                "posterior_mean": float(
                    np.mean(
                        values
                    )
                ),
                "posterior_median": float(
                    np.median(
                        values
                    )
                ),
                "ci50_lower": float(
                    np.quantile(
                        values,
                        0.25,
                    )
                ),
                "ci50_upper": float(
                    np.quantile(
                        values,
                        0.75,
                    )
                ),
                "ci95_lower": float(
                    np.quantile(
                        values,
                        0.025,
                    )
                ),
                "ci95_upper": float(
                    np.quantile(
                        values,
                        0.975,
                    )
                ),
                "probability_rank_1": float(
                    np.mean(
                        rank_matrix[
                            :,
                            mode_index
                        ] == 1
                    )
                ),
                "probability_rank_top_2": float(
                    np.mean(
                        rank_matrix[
                            :,
                            mode_index
                        ] <= 2
                    )
                ),
                "mean_rank": float(
                    np.mean(
                        rank_matrix[
                            :,
                            mode_index
                        ]
                    )
                ),
                "n_observed_supporting_structures": (
                    observed_supporting_structures
                ),
            }
        )

    return (
        pd.DataFrame(
            rows
        )
        .sort_values(
            "posterior_median",
            ascending=False,
        )
        .reset_index(
            drop=True
        )
    )


def dominance_probability_matrix(
    draws: np.ndarray,
    mode_labels: list[str],
) -> pd.DataFrame:
    n_modes = len(
        mode_labels
    )

    matrix = np.zeros(
        (
            n_modes,
            n_modes,
        ),
        dtype=float,
    )

    for first_index in range(
        n_modes
    ):
        for second_index in range(
            n_modes
        ):
            if first_index == second_index:
                matrix[
                    first_index,
                    second_index,
                ] = 0.5

            else:
                matrix[
                    first_index,
                    second_index,
                ] = float(
                    np.mean(
                        draws[
                            :,
                            first_index
                        ]
                        > draws[
                            :,
                            second_index
                        ]
                    )
                )

    return pd.DataFrame(
        matrix,
        index=mode_labels,
        columns=mode_labels,
    )


def rank_probability_table(
    draws: np.ndarray,
    mode_labels: list[str],
) -> pd.DataFrame:
    n_draws, n_modes = draws.shape

    order_per_draw = np.argsort(
        -draws,
        axis=1,
    )

    rows = []

    for rank_index in range(
        n_modes
    ):
        counts = np.bincount(
            order_per_draw[
                :,
                rank_index
            ],
            minlength=n_modes,
        )

        for mode_index, count in enumerate(
            counts
        ):
            rows.append(
                {
                    "mode": mode_labels[
                        mode_index
                    ],
                    "rank": rank_index + 1,
                    "probability": float(
                        count
                        / n_draws
                    ),
                }
            )

    return pd.DataFrame(
        rows
    )


# ============================================================
# Plotting
# ============================================================

def _select_display_modes(
    summary: pd.DataFrame,
    top_n: int | None,
    minimum_median_support: float,
) -> list[str]:
    selected = summary[
        summary[
            "posterior_median"
        ]
        >= minimum_median_support
    ][
        "mode"
    ].tolist()

    if top_n is not None:
        selected = selected[
            :top_n
        ]

    if not selected:
        selected = summary[
            "mode"
        ].head(
            top_n or len(
                summary
            )
        ).tolist()

    return selected


def save_support_violin_plot(
    equal_draws: np.ndarray,
    adjusted_draws: np.ndarray,
    mode_labels: list[str],
    summary: pd.DataFrame,
    output_file: Path,
    title: str,
    top_n: int | None,
    minimum_median_support: float,
):
    selected_modes = _select_display_modes(
        summary=summary,
        top_n=top_n,
        minimum_median_support=(
            minimum_median_support
        ),
    )

    mode_to_index = {
        mode_label: index
        for index, mode_label in enumerate(
            mode_labels
        )
    }

    selected_modes = [
        mode
        for mode in selected_modes
        if mode in mode_to_index
    ]

    figure_height = max(
        5.0,
        0.48
        * len(
            selected_modes
        )
        + 2.0,
    )

    figure, axis = plt.subplots(
        figsize=(
            10.0,
            figure_height,
        )
    )

    positions_equal = np.arange(
        len(
            selected_modes
        )
    ) * 2.0

    positions_adjusted = (
        positions_equal
        + 0.65
    )

    equal_values = [
        equal_draws[
            :,
            mode_to_index[
                mode
            ],
        ]
        for mode in selected_modes
    ]

    adjusted_values = [
        adjusted_draws[
            :,
            mode_to_index[
                mode
            ],
        ]
        for mode in selected_modes
    ]

    equal_violin = axis.violinplot(
        equal_values,
        positions=positions_equal,
        vert=False,
        widths=0.55,
        showmeans=False,
        showmedians=True,
        showextrema=False,
    )

    adjusted_violin = axis.violinplot(
        adjusted_values,
        positions=positions_adjusted,
        vert=False,
        widths=0.55,
        showmeans=False,
        showmedians=True,
        showextrema=False,
    )

    for body in equal_violin[
        "bodies"
    ]:
        body.set_facecolor(
            "#7aa6c2"
        )

        body.set_alpha(
            0.70
        )

    for body in adjusted_violin[
        "bodies"
    ]:
        body.set_facecolor(
            "#d98c4c"
        )

        body.set_alpha(
            0.70
        )

    axis.set_yticks(
        positions_equal
        + 0.325
    )

    axis.set_yticklabels(
        selected_modes
    )

    axis.set_xlabel(
        "Posterior ensemble support"
    )

    axis.set_ylabel(
        "Mode"
    )

    axis.set_xlim(
        0.0,
        1.0,
    )

    axis.set_title(
        title
    )

    axis.scatter(
        [],
        [],
        color="#7aa6c2",
        label="Equal structure weighting",
    )

    axis.scatter(
        [],
        [],
        color="#d98c4c",
        label="Similarity-adjusted weighting",
    )

    axis.legend(
        loc="best"
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


def save_dominance_heatmap(
    dominance_df: pd.DataFrame,
    output_file: Path,
    title: str,
    maximum_modes: int = 20,
):
    if len(
        dominance_df
    ) > maximum_modes:
        selected = dominance_df.index[
            :maximum_modes
        ]

        dominance_df = dominance_df.loc[
            selected,
            selected,
        ]

    figure_size = max(
        7.0,
        0.52
        * len(
            dominance_df
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
        dominance_df.to_numpy(
            dtype=float
        ),
        vmin=0.0,
        vmax=1.0,
        cmap="coolwarm",
        aspect="auto",
    )

    axis.set_xticks(
        np.arange(
            len(
                dominance_df.columns
            )
        )
    )

    axis.set_xticklabels(
        dominance_df.columns,
        rotation=90,
    )

    axis.set_yticks(
        np.arange(
            len(
                dominance_df.index
            )
        )
    )

    axis.set_yticklabels(
        dominance_df.index
    )

    axis.set_xlabel(
        "Comparison mode"
    )

    axis.set_ylabel(
        "Row mode"
    )

    axis.set_title(
        title
    )

    colour_bar = figure.colorbar(
        image,
        ax=axis,
    )

    colour_bar.set_label(
        "P(row mode support > column mode support)"
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


def save_rank_probability_plot(
    rank_df: pd.DataFrame,
    summary: pd.DataFrame,
    output_file: Path,
    title: str,
    top_n_modes: int = 15,
    maximum_rank: int = 5,
):
    selected_modes = summary[
        "mode"
    ].head(
        top_n_modes
    ).tolist()

    filtered = rank_df[
        rank_df[
            "mode"
        ].isin(
            selected_modes
        )
        & (
            rank_df[
                "rank"
            ]
            <= maximum_rank
        )
    ]

    matrix = (
        filtered
        .pivot(
            index="mode",
            columns="rank",
            values="probability",
        )
        .fillna(
            0.0
        )
        .reindex(
            selected_modes
        )
    )

    figure_height = max(
        5.0,
        0.45
        * len(
            matrix
        )
        + 2.0,
    )

    figure, axis = plt.subplots(
        figsize=(
            10.0,
            figure_height,
        )
    )

    left = np.zeros(
        len(
            matrix
        ),
        dtype=float,
    )

    for rank in matrix.columns:
        values = matrix[
            rank
        ].to_numpy(
            dtype=float
        )

        axis.barh(
            np.arange(
                len(
                    matrix
                )
            ),
            values,
            left=left,
            label=f"Rank {int(rank)}",
        )

        left += values

    axis.set_yticks(
        np.arange(
            len(
                matrix
            )
        )
    )

    axis.set_yticklabels(
        matrix.index
    )

    axis.set_xlim(
        0.0,
        1.0,
    )

    axis.set_xlabel(f"Posterior probability of rank 1-{maximum_rank}"
    )

    axis.set_title(
        title
    )

    axis.legend(
        loc="best"
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


def save_structure_weight_plot(
    weight_draws: np.ndarray,
    source_names: list[str],
    output_file: Path,
):
    figure_height = max(
        5.0,
        0.55
        * len(
            source_names
        )
        + 2.0,
    )

    figure, axis = plt.subplots(
        figsize=(
            9.0,
            figure_height,
        )
    )

    values = [
        weight_draws[
            :,
            structure_index
        ]
        for structure_index in range(
            len(
                source_names
            )
        )
    ]

    positions = np.arange(
        len(
            source_names
        )
    )

    violins = axis.violinplot(
        values,
        positions=positions,
        vert=False,
        widths=0.75,
        showmeans=False,
        showmedians=True,
        showextrema=False,
    )

    for body in violins[
        "bodies"
    ]:
        body.set_facecolor(
            "#8c6bb1"
        )

        body.set_alpha(
            0.75
        )

    axis.set_yticks(
        positions
    )

    axis.set_yticklabels(
        source_names
    )

    axis.set_xlim(
        0.0,
        max(
            0.5,
            float(
                np.max(
                    weight_draws
                )
            )
            * 1.05,
        ),
    )

    axis.set_xlabel(
        "Learned campaign weight"
    )

    axis.set_ylabel(
        "Structure-level GlideEM ensemble"
    )

    axis.set_title(
        "Posterior distribution of similarity-adjusted "
        "structure weights"
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


def save_js_heatmap(
    matrix: np.ndarray,
    source_names: list[str],
    output_file: Path,
):
    figure_size = max(
        7.0,
        0.60
        * len(
            source_names
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
        matrix,
        vmin=0.0,
        vmax=max(
            1.0,
            float(
                np.max(
                    matrix
                )
            ),
        ),
        cmap="magma",
    )

    axis.set_xticks(
        np.arange(
            len(
                source_names
            )
        )
    )

    axis.set_xticklabels(
        source_names,
        rotation=90,
    )

    axis.set_yticks(
        np.arange(
            len(
                source_names
            )
        )
    )

    axis.set_yticklabels(
        source_names
    )

    axis.set_title(
        "Posterior mean structure-landscape "
        "Jensen-Shannon distance"
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


# ============================================================
# Main entry point
# ============================================================

def run_ensemble_support_analysis(
    pose_table: pd.DataFrame,
    output_dir: Path,
    config: dict,
) -> dict:
    _validate_pose_table(
        pose_table
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

    n_draws = int(
        config.get(
            "n_draws",
            10000,
        )
    )

    random_seed = int(
        config.get(
            "random_seed",
            20261023,
        )
    )

    alpha = float(
        config.get(
            "dirichlet_alpha",
            0.5,
        )
    )

    support_threshold = float(
        config.get(
            "support_threshold",
            0.10,
        )
    )

    effective_cfg = config.get(
        "effective_counts",
        {},
    )

    effective_enabled = bool(
        effective_cfg.get(
            "enabled",
            True,
        )
    )

    effective_minimum = float(
        effective_cfg.get(
            "minimum",
            1.0,
        )
    )

    redundancy_cfg = config.get(
        "redundancy_adjustment",
        {},
    )

    bandwidth_multipliers = [
        float(
            value
        )
        for value in redundancy_cfg.get(
            "sensitivity_multipliers",
            [
                0.5,
                1.0,
                2.0,
            ],
        )
    ]

    if 1.0 not in bandwidth_multipliers:
        bandwidth_multipliers.append(
            1.0
        )

        bandwidth_multipliers = sorted(
            bandwidth_multipliers
        )

    display_cfg = config.get(
        "display",
        {},
    )

    top_n_subfamilies = int(
        display_cfg.get(
            "top_n_subfamilies",
            20,
        )
    )

    minimum_median_support = float(
        display_cfg.get(
            "minimum_median_support",
            0.02,
        )
    )

    if n_draws < 100:
        raise ValueError(
            "ensemble_support.n_draws must be at least 100."
        )

    if alpha <= 0.0:
        raise ValueError(
            "ensemble_support.dirichlet_alpha must be positive."
        )

    rng = np.random.default_rng(
        random_seed
    )

    source_names = sorted(
        pose_table[
            "source_name"
        ].unique()
        .tolist()
    )

    spatial_labels = _mode_labels(
        pose_table,
        mode_level="spatial",
    )

    joint_labels = _mode_labels(
        pose_table,
        mode_level="joint",
    )

    effective_sizes_df = structure_effective_sizes(
        pose_table=pose_table,
        enabled=effective_enabled,
        minimum=effective_minimum,
    )

    effective_sizes_df = (
        effective_sizes_df
        .set_index(
            "source_name"
        )
        .loc[
            source_names
        ]
        .reset_index()
    )

    effective_sizes = effective_sizes_df[
        "effective_nposes"
    ].to_numpy(
        dtype=float
    )

    effective_sizes_df.to_csv(
        output_dir
        / "structure_effective_sample_sizes.csv",
        index=False,
    )

    spatial_observed = observed_mode_proportions(
        pose_table=pose_table,
        mode_level="spatial",
        source_names=source_names,
        mode_labels=spatial_labels,
    )

    joint_observed = observed_mode_proportions(
        pose_table=pose_table,
        mode_level="joint",
        source_names=source_names,
        mode_labels=joint_labels,
    )

    spatial_structure_draws = draw_structure_posteriors(
        observed_proportions=spatial_observed,
        effective_sizes=effective_sizes,
        n_draws=n_draws,
        alpha=alpha,
        rng=rng,
    )

    joint_structure_draws = draw_structure_posteriors(
        observed_proportions=joint_observed,
        effective_sizes=effective_sizes,
        n_draws=n_draws,
        alpha=alpha,
        rng=rng,
    )

    # Learn similarity entirely from the richer joint landscape.
    joint_campaign = campaign_support_draws(
        structure_draws=joint_structure_draws,
        bandwidth_multipliers=(
            bandwidth_multipliers
        ),
    )

    spatial_equal = np.mean(
        spatial_structure_draws,
        axis=1,
    ).astype(
        np.float32
    )

    spatial_adjusted = apply_learned_weights(
        structure_draws=spatial_structure_draws,
        weights_by_multiplier=(
            joint_campaign[
                "weights"
            ]
        ),
    )

    primary_multiplier = 1.0

    joint_equal = joint_campaign[
        "equal"
    ]

    joint_primary = joint_campaign[
        "adjusted"
    ][
        primary_multiplier
    ]

    spatial_primary = spatial_adjusted[
        primary_multiplier
    ]

    # --------------------------------------------------
    # Save posterior draws
    # --------------------------------------------------

    draw_payload = {
        "source_names": np.asarray(
            source_names,
            dtype=str,
        ),
        "spatial_labels": np.asarray(
            spatial_labels,
            dtype=str,
        ),
        "joint_labels": np.asarray(
            joint_labels,
            dtype=str,
        ),
        "spatial_structure_draws": (
            spatial_structure_draws
        ),
        "joint_structure_draws": (
            joint_structure_draws
        ),
        "spatial_equal_draws": (
            spatial_equal
        ),
        "joint_equal_draws": (
            joint_equal
        ),
        "mean_js_matrix": (
            joint_campaign[
                "mean_js_matrix"
            ]
        ),
    }

    for multiplier in bandwidth_multipliers:
        key_suffix = str(
            multiplier
        ).replace(
            ".",
            "p",
        )

        draw_payload[
            f"spatial_adjusted_{key_suffix}"
        ] = spatial_adjusted[
            multiplier
        ]

        draw_payload[
            f"joint_adjusted_{key_suffix}"
        ] = joint_campaign[
            "adjusted"
        ][
            multiplier
        ]

        draw_payload[
            f"structure_weights_{key_suffix}"
        ] = joint_campaign[
            "weights"
        ][
            multiplier
        ]

        draw_payload[
            f"bandwidths_{key_suffix}"
        ] = joint_campaign[
            "bandwidths"
        ][
            multiplier
        ]

    np.savez_compressed(
        output_dir
        / "ensemble_support_posterior_draws.npz",
        **draw_payload,
    )

    # --------------------------------------------------
    # Summaries
    # --------------------------------------------------

    spatial_equal_summary = posterior_summary(
        draws=spatial_equal,
        mode_labels=spatial_labels,
        observed_proportions=spatial_observed,
        support_threshold=support_threshold,
        analysis_name="spatial_equal",
    )

    spatial_adjusted_summary = posterior_summary(
        draws=spatial_primary,
        mode_labels=spatial_labels,
        observed_proportions=spatial_observed,
        support_threshold=support_threshold,
        analysis_name="spatial_similarity_adjusted",
    )

    joint_equal_summary = posterior_summary(
        draws=joint_equal,
        mode_labels=joint_labels,
        observed_proportions=joint_observed,
        support_threshold=support_threshold,
        analysis_name="joint_equal",
    )

    joint_adjusted_summary = posterior_summary(
        draws=joint_primary,
        mode_labels=joint_labels,
        observed_proportions=joint_observed,
        support_threshold=support_threshold,
        analysis_name="joint_similarity_adjusted",
    )

    # hierarchical_summary = (
    #     build_hierarchical_support_summary(
    #         spatial_adjusted_summary=(
    #             spatial_adjusted_summary
    #         ),
    #         conditional_results=(
    #             conditional_results
    #         ),
    #         joint_adjusted_summary=(
    #             joint_adjusted_summary
    #         ),
    #     )
    # )

    # hierarchical_summary.to_csv(
    #     output_dir
    #     / "hierarchical_mode_support.csv",
    #     index=False,
    # )

    # save_conditional_support_plot(
    #     conditional_results=(
    #         conditional_results
    #     ),
    #     output_file=(
    #         output_dir
    #         / "conditional_subfamily_support.png"
    #     ),
    # )

    spatial_summary = pd.concat(
        [
            spatial_equal_summary,
            spatial_adjusted_summary,
        ],
        ignore_index=True,
    )

    joint_summary = pd.concat(
        [
            joint_equal_summary,
            joint_adjusted_summary,
        ],
        ignore_index=True,
    )

    spatial_summary.to_csv(
        output_dir
        / "spatial_support_posterior.csv",
        index=False,
    )

    joint_summary.to_csv(
        output_dir
        / "subfamily_support_posterior.csv",
        index=False,
    )

    # --------------------------------------------------
    # Weight and bandwidth summaries
    # --------------------------------------------------

    primary_weights = joint_campaign[
        "weights"
    ][
        primary_multiplier
    ]

    weight_rows = []

    for structure_index, source_name in enumerate(
        source_names
    ):
        values = primary_weights[
            :,
            structure_index
        ]

        weight_rows.append(
            {
                "source_name": source_name,
                "posterior_mean_weight": float(
                    np.mean(
                        values
                    )
                ),
                "posterior_median_weight": float(
                    np.median(
                        values
                    )
                ),
                "ci50_lower": float(
                    np.quantile(
                        values,
                        0.25,
                    )
                ),
                "ci50_upper": float(
                    np.quantile(
                        values,
                        0.75,
                    )
                ),
                "ci95_lower": float(
                    np.quantile(
                        values,
                        0.025,
                    )
                ),
                "ci95_upper": float(
                    np.quantile(
                        values,
                        0.975,
                    )
                ),
            }
        )

    pd.DataFrame(
        weight_rows
    ).to_csv(
        output_dir
        / "structure_similarity_weights.csv",
        index=False,
    )

    # --------------------------------------------------
    # Conditional orientation/contact support
    #
    # Calculate:
    #
    #     P(subfamily | parent spatial family)
    #
    # This must occur after primary_weights has been created.
    # --------------------------------------------------

    conditional_results = (
        run_conditional_subfamily_analysis(
            pose_table=pose_table,
            source_names=source_names,
            effective_sizes=effective_sizes,
            spatial_labels=spatial_labels,
            spatial_structure_draws=(
                spatial_structure_draws
            ),
            campaign_structure_weights=(
                primary_weights
            ),
            n_draws=n_draws,
            alpha=alpha,
            support_threshold=(
                support_threshold
            ),
            rng=rng,
        )
    )

    if not conditional_results:
        raise RuntimeError(
            "Conditional subfamily analysis produced no results."
        )

    expected_conditional_families = set(
        spatial_labels
    )

    observed_conditional_families = set(
        conditional_results
    )

    if (
        observed_conditional_families
        != expected_conditional_families
    ):
        raise RuntimeError(
            "Conditional family labels do not match spatial "
            "posterior labels. "
            f"Expected: {sorted(expected_conditional_families)}; "
            f"observed: {sorted(observed_conditional_families)}"
        )    

    # --------------------------------------------------
    # Save per-family and combined conditional summaries
    # --------------------------------------------------

    conditional_summary_tables = []
    conditional_draw_payload = {}

    for family_label, result in sorted(
        conditional_results.items(),
        key=lambda item: int(
            item[
                0
            ][
                1:
            ]
        ),
    ):
        family_safe = (
            family_label.lower()
        )

        family_summary = pd.concat(
            [
                result[
                    "summary_equal"
                ],
                result[
                    "summary_adjusted"
                ],
            ],
            ignore_index=True,
        )

        family_summary.to_csv(
            output_dir
            / (
                f"{family_safe}"
                "_conditional_subfamily_support.csv"
            ),
            index=False,
        )

        conditional_summary_tables.append(
            family_summary
        )

        conditional_draw_payload[
            f"{family_safe}_labels"
        ] = np.asarray(
            result[
                "subfamily_labels"
            ],
            dtype=str,
        )

        conditional_draw_payload[
            f"{family_safe}_equal_draws"
        ] = result[
            "equal_draws"
        ]

        conditional_draw_payload[
            f"{family_safe}_adjusted_draws"
        ] = result[
            "adjusted_draws"
        ]

        conditional_draw_payload[
            f"{family_safe}_structure_draws"
        ] = result[
            "structure_draws"
        ]

        conditional_draw_payload[
            f"{family_safe}_conditional_weights"
        ] = result[
            "conditional_weights"
        ]

    if not conditional_summary_tables:
        raise RuntimeError(
            "No conditional subfamily summaries were generated."
        )

    all_conditional_summary = pd.concat(
        conditional_summary_tables,
        ignore_index=True,
    )

    all_conditional_summary.to_csv(
        output_dir
        / "conditional_subfamily_support.csv",
        index=False,
    )

    np.savez_compressed(
        output_dir
        / "conditional_subfamily_posterior_draws.npz",
        **conditional_draw_payload,
    )

    # --------------------------------------------------
    # Hierarchical summary:
    #
    #     P(spatial family)
    #     P(subfamily | spatial family)
    #     P(spatial family and subfamily)
    # --------------------------------------------------

    hierarchical_summary = (
        build_hierarchical_support_summary(
            spatial_adjusted_summary=(
                spatial_adjusted_summary
            ),
            conditional_results=(
                conditional_results
            ),
            joint_adjusted_summary=(
                joint_adjusted_summary
            ),
        )
    )

    hierarchical_summary.to_csv(
        output_dir
        / "hierarchical_mode_support.csv",
        index=False,
    )

    # --------------------------------------------------
    # Conditional posterior plot
    # --------------------------------------------------

    save_conditional_support_plot(
        conditional_results=(
            conditional_results
        ),
        output_file=(
            output_dir
            / "conditional_subfamily_support.png"
        ),
    )

    logger.info(
        "Generated conditional subfamily posterior analysis "
        "for %d spatial families",
        len(
            conditional_results
        ),
    )

    sensitivity_rows = []

    for multiplier in bandwidth_multipliers:
        bandwidth_values = joint_campaign[
            "bandwidths"
        ][
            multiplier
        ]

        sensitivity_rows.append(
            {
                "bandwidth_multiplier": multiplier,
                "median_bandwidth": float(
                    np.median(
                        bandwidth_values
                    )
                ),
                "ci95_lower": float(
                    np.quantile(
                        bandwidth_values,
                        0.025,
                    )
                ),
                "ci95_upper": float(
                    np.quantile(
                        bandwidth_values,
                        0.975,
                    )
                ),
            }
        )

    pd.DataFrame(
        sensitivity_rows
    ).to_csv(
        output_dir
        / "similarity_bandwidth_sensitivity.csv",
        index=False,
    )

    # --------------------------------------------------
    # Dominance and ranks
    # --------------------------------------------------

    spatial_dominance = dominance_probability_matrix(
        draws=spatial_primary,
        mode_labels=spatial_labels,
    )

    joint_dominance = dominance_probability_matrix(
        draws=joint_primary,
        mode_labels=joint_labels,
    )

    spatial_dominance.to_csv(
        output_dir
        / "spatial_dominance_probabilities.csv"
    )

    joint_dominance.to_csv(
        output_dir
        / "subfamily_dominance_probabilities.csv"
    )

    spatial_ranks = rank_probability_table(
        draws=spatial_primary,
        mode_labels=spatial_labels,
    )

    joint_ranks = rank_probability_table(
        draws=joint_primary,
        mode_labels=joint_labels,
    )

    spatial_ranks.to_csv(
        output_dir
        / "spatial_rank_probabilities.csv",
        index=False,
    )

    joint_ranks.to_csv(
        output_dir
        / "subfamily_rank_probabilities.csv",
        index=False,
    )

    # --------------------------------------------------
    # Plots
    # --------------------------------------------------

    save_support_violin_plot(
        equal_draws=spatial_equal,
        adjusted_draws=spatial_primary,
        mode_labels=spatial_labels,
        summary=spatial_adjusted_summary,
        output_file=(
            output_dir
            / "spatial_support_distributions.png"
        ),
        title=(
            "Posterior ensemble support: spatial families"
        ),
        top_n=None,
        minimum_median_support=0.0,
    )

    save_support_violin_plot(
        equal_draws=joint_equal,
        adjusted_draws=joint_primary,
        mode_labels=joint_labels,
        summary=joint_adjusted_summary,
        output_file=(
            output_dir
            / "subfamily_support_distributions.png"
        ),
        title=(
            "Posterior ensemble support: "
            "orientation/contact subfamilies"
        ),
        top_n=top_n_subfamilies,
        minimum_median_support=(
            minimum_median_support
        ),
    )

    save_dominance_heatmap(
        dominance_df=spatial_dominance,
        output_file=(
            output_dir
            / "spatial_dominance_heatmap.png"
        ),
        title=(
            "Posterior dominance probabilities: "
            "spatial families"
        ),
        maximum_modes=len(
            spatial_labels
        ),
    )

    selected_joint_modes = (
        joint_adjusted_summary[
            "mode"
        ]
        .head(
            top_n_subfamilies
        )
        .tolist()
    )

    joint_dominance_display = (
        joint_dominance.loc[
            selected_joint_modes,
            selected_joint_modes,
        ]
    )

    save_dominance_heatmap(
        dominance_df=joint_dominance_display,
        output_file=(
            output_dir
            / "subfamily_dominance_heatmap.png"
        ),
        title=(
            "Posterior dominance probabilities: "
            "top orientation/contact subfamilies"
        ),
        maximum_modes=top_n_subfamilies,
    )

    save_rank_probability_plot(
        rank_df=spatial_ranks,
        summary=spatial_adjusted_summary,
        output_file=(
            output_dir
            / "spatial_rank_probabilities.png"
        ),
        title=(
            "Posterior rank probabilities: spatial families"
        ),
        top_n_modes=len(
            spatial_labels
        ),
        maximum_rank=min(
            5,
            len(
                spatial_labels
            ),
        ),
    )

    save_rank_probability_plot(
        rank_df=joint_ranks,
        summary=joint_adjusted_summary,
        output_file=(
            output_dir
            / "subfamily_rank_probabilities.png"
        ),
        title=(
            "Posterior rank probabilities: "
            "orientation/contact subfamilies"
        ),
        top_n_modes=top_n_subfamilies,
        maximum_rank=5,
    )

    save_structure_weight_plot(
        weight_draws=primary_weights,
        source_names=source_names,
        output_file=(
            output_dir
            / "structure_weight_distributions.png"
        ),
    )

    save_js_heatmap(
        matrix=joint_campaign[
            "mean_js_matrix"
        ],
        source_names=source_names,
        output_file=(
            output_dir
            / "posterior_mean_structure_js_heatmap.png"
        ),
    )

    metadata = {
        "n_draws": n_draws,
        "random_seed": random_seed,
        "dirichlet_alpha": alpha,
        "effective_counts_enabled": (
            effective_enabled
        ),
        "support_threshold": (
            support_threshold
        ),
        "primary_bandwidth_multiplier": 1.0,
        "bandwidth_multipliers": (
            bandwidth_multipliers
        ),
        "n_structures": len(
            source_names
        ),
        "n_spatial_families": len(
            spatial_labels
        ),
        "n_joint_subfamilies": len(
            joint_labels
        ),
    }

    (
        output_dir
        / "ensemble_support_metadata.json"
    ).write_text(
        json.dumps(
            metadata,
            indent=2,
        )
    )

    return {
        "posterior_draws": (
            output_dir
            / "ensemble_support_posterior_draws.npz"
        ),
        "spatial_summary": (
            output_dir
            / "spatial_support_posterior.csv"
        ),
        "subfamily_summary": (
            output_dir
            / "subfamily_support_posterior.csv"
        ),
        "structure_weights": (
            output_dir
            / "structure_similarity_weights.csv"
        ),
        "spatial_distributions_plot": (
            output_dir
            / "spatial_support_distributions.png"
        ),
        "subfamily_distributions_plot": (
            output_dir
            / "subfamily_support_distributions.png"
        ),
        "conditional_subfamily_summary": (
            output_dir
            / "conditional_subfamily_support.csv"
        ),
        "hierarchical_mode_support": (
            output_dir
            / "hierarchical_mode_support.csv"
        ),
        "conditional_subfamily_plot": (
            output_dir
            / "conditional_subfamily_support.png"
        ),
        "conditional_posterior_draws": (
            output_dir
            / "conditional_subfamily_posterior_draws.npz"
        ),
        "_joint_adjusted_summary_df": (
            joint_adjusted_summary
        ),
    }

def draw_conditional_subfamily_posteriors(
    pose_table: pd.DataFrame,
    source_names: list[str],
    effective_sizes: np.ndarray,
    n_draws: int,
    alpha: float,
    rng: np.random.Generator,
) -> dict:
    """
    Draw P(subfamily | spatial family, structure).

    Returns a dictionary keyed by spatial-family label, for example:

        {
            "F1": {
                "subfamily_labels": ["F1.1", "F1.2", ...],
                "structure_draws": array(
                    shape=(
                        n_draws,
                        n_structures,
                        n_subfamilies,
                    )
                ),
                "observed_family_fraction": array(
                    shape=(n_structures,)
                ),
                "observed_conditional_proportions": array(
                    shape=(
                        n_structures,
                        n_subfamilies,
                    )
                ),
            },
            ...
        }

    The effective sampling depth available to a family is scaled
    by the observed fraction of the structure's poses occupying
    that spatial family.
    """
    working = pose_table.copy()

    working[
        "_spatial_family_id"
    ] = working[
        "global_spatial_family_id"
    ].astype(
        int
    )

    working[
        "_subfamily_id"
    ] = working[
        "orientation_contact_subfamily_id"
    ].astype(
        int
    )

    source_to_index = {
        source_name: index
        for index, source_name in enumerate(
            source_names
        )
    }

    structure_totals = (
        working
        .groupby(
            "source_name"
        )
        .size()
        .to_dict()
    )

    result = {}

    spatial_family_ids = sorted(
        working[
            "_spatial_family_id"
        ].unique()
        .tolist()
    )

    for family_id in spatial_family_ids:
        family_label = (
            f"F{int(family_id)}"
        )

        family_data = working[
            working[
                "_spatial_family_id"
            ] == family_id
        ]

        subfamily_ids = sorted(
            family_data[
                "_subfamily_id"
            ].unique()
            .tolist()
        )

        subfamily_labels = [
            (
                f"F{int(family_id)}."
                f"{int(subfamily_id)}"
            )
            for subfamily_id in subfamily_ids
        ]

        subfamily_to_index = {
            int(
                subfamily_id
            ): index
            for index, subfamily_id in enumerate(
                subfamily_ids
            )
        }

        n_structures = len(
            source_names
        )

        n_subfamilies = len(
            subfamily_ids
        )

        observed_family_fraction = np.zeros(
            n_structures,
            dtype=float,
        )

        observed_conditional = np.zeros(
            (
                n_structures,
                n_subfamilies,
            ),
            dtype=float,
        )

        structure_draws = np.zeros(
            (
                n_draws,
                n_structures,
                n_subfamilies,
            ),
            dtype=np.float32,
        )

        for source_name in source_names:
            structure_index = source_to_index[
                source_name
            ]

            source_total = int(
                structure_totals.get(
                    source_name,
                    0,
                )
            )

            source_family_data = family_data[
                family_data[
                    "source_name"
                ] == source_name
            ]

            family_count = int(
                len(
                    source_family_data
                )
            )

            if source_total > 0:
                family_fraction = float(
                    family_count
                    / source_total
                )

            else:
                family_fraction = 0.0

            observed_family_fraction[
                structure_index
            ] = family_fraction

            counts = np.zeros(
                n_subfamilies,
                dtype=float,
            )

            if family_count > 0:
                subfamily_counts = (
                    source_family_data[
                        "_subfamily_id"
                    ]
                    .value_counts()
                )

                for (
                    subfamily_id,
                    count,
                ) in subfamily_counts.items():
                    counts[
                        subfamily_to_index[
                            int(
                                subfamily_id
                            )
                        ]
                    ] = float(
                        count
                    )

                observed_conditional[
                    structure_index
                ] = (
                    counts
                    / float(
                        family_count
                    )
                )

            # Allocate this structure's effective sampling depth
            # to the family in proportion to the amount of that
            # structure's ensemble occupying the family.
            family_effective_n = (
                float(
                    effective_sizes[
                        structure_index
                    ]
                )
                * family_fraction
            )

            # No evidence for this spatial family from this
            # structure: retain only the symmetric prior.
            if family_count == 0:
                effective_counts = np.zeros(
                    n_subfamilies,
                    dtype=float,
                )

            else:
                effective_counts = (
                    observed_conditional[
                        structure_index
                    ]
                    * family_effective_n
                )

            parameters = (
                effective_counts
                + alpha
            )

            structure_draws[
                :,
                structure_index,
                :,
            ] = rng.dirichlet(
                parameters,
                size=n_draws,
            ).astype(
                np.float32
            )

        result[
            family_label
        ] = {
            "family_id": int(
                family_id
            ),
            "subfamily_ids": [
                int(
                    value
                )
                for value in subfamily_ids
            ],
            "subfamily_labels": (
                subfamily_labels
            ),
            "structure_draws": (
                structure_draws
            ),
            "observed_family_fraction": (
                observed_family_fraction
            ),
            "observed_conditional_proportions": (
                observed_conditional
            ),
        }

    return result

def conditional_campaign_draws(
    conditional_structure_draws: np.ndarray,
    parent_spatial_structure_draws: np.ndarray,
    campaign_structure_weights: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Combine structure-level conditional subfamily posteriors.

    Structures are weighted by both:

        1. objectively learned campaign redundancy weight;
        2. that draw's support for the parent spatial family.

    This prevents a structure that never samples F1 from voting
    strongly on P(F1.x | F1).

    Parameters
    ----------
    conditional_structure_draws
        Shape:
            n_draws x n_structures x n_subfamilies

    parent_spatial_structure_draws
        Shape:
            n_draws x n_structures

        Posterior P(parent spatial family | structure).

    campaign_structure_weights
        Shape:
            n_draws x n_structures

    Returns
    -------
    campaign_draws
        n_draws x n_subfamilies

    conditional_weights
        n_draws x n_structures
    """
    raw_conditional_weights = (
        campaign_structure_weights
        * parent_spatial_structure_draws
    )

    row_sums = np.sum(
        raw_conditional_weights,
        axis=1,
        keepdims=True,
    )

    conditional_weights = np.divide(
        raw_conditional_weights,
        row_sums,
        out=np.zeros_like(
            raw_conditional_weights,
            dtype=float,
        ),
        where=row_sums > 0.0,
    )

    zero_rows = (
        row_sums[
            :,
            0
        ] <= 0.0
    )

    if np.any(
        zero_rows
    ):
        conditional_weights[
            zero_rows
        ] = (
            1.0
            / conditional_weights.shape[
                1
            ]
        )

    campaign_draws = np.einsum(
        "bs,bsm->bm",
        conditional_weights,
        conditional_structure_draws,
    ).astype(
        np.float32
    )

    return (
        campaign_draws,
        conditional_weights.astype(
            np.float32
        ),
    )

def conditional_equal_structure_draws(
    conditional_structure_draws: np.ndarray,
    parent_spatial_structure_draws: np.ndarray,
) -> np.ndarray:
    """
    Conditional campaign support without landscape-similarity
    adjustment.

    Structures are still weighted by posterior support for the
    parent spatial family.
    """
    weights = np.asarray(
        parent_spatial_structure_draws,
        dtype=float,
    )

    row_sums = weights.sum(
        axis=1,
        keepdims=True,
    )

    weights = np.divide(
        weights,
        row_sums,
        out=np.zeros_like(
            weights,
            dtype=float,
        ),
        where=row_sums > 0.0,
    )

    zero_rows = (
        row_sums[
            :,
            0
        ] <= 0.0
    )

    if np.any(
        zero_rows
    ):
        weights[
            zero_rows
        ] = (
            1.0
            / weights.shape[
                1
            ]
        )

    return np.einsum(
        "bs,bsm->bm",
        weights,
        conditional_structure_draws,
    ).astype(
        np.float32
    )

def run_conditional_subfamily_analysis(
    pose_table: pd.DataFrame,
    source_names: list[str],
    effective_sizes: np.ndarray,
    spatial_labels: list[str],
    spatial_structure_draws: np.ndarray,
    campaign_structure_weights: np.ndarray,
    n_draws: int,
    alpha: float,
    support_threshold: float,
    rng: np.random.Generator,
) -> dict:
    """
    Calculate P(subfamily | spatial family) for every spatial family.
    """
    raw_results = draw_conditional_subfamily_posteriors(
        pose_table=pose_table,
        source_names=source_names,
        effective_sizes=effective_sizes,
        n_draws=n_draws,
        alpha=alpha,
        rng=rng,
    )

    spatial_label_to_index = {
        label: index
        for index, label in enumerate(
            spatial_labels
        )
    }

    output = {}

    for family_label, family_result in (
        raw_results.items()
    ):
        if family_label not in spatial_label_to_index:
            raise RuntimeError(
                f"Conditional family {family_label} does not "
                "exist in spatial posterior labels."
            )

        parent_index = spatial_label_to_index[
            family_label
        ]

        parent_structure_draws = (
            spatial_structure_draws[
                :,
                :,
                parent_index,
            ]
        )

        conditional_draws = family_result[
            "structure_draws"
        ]

        equal_draws = (
            conditional_equal_structure_draws(
                conditional_structure_draws=(
                    conditional_draws
                ),
                parent_spatial_structure_draws=(
                    parent_structure_draws
                ),
            )
        )

        (
            adjusted_draws,
            conditional_weights,
        ) = conditional_campaign_draws(
            conditional_structure_draws=(
                conditional_draws
            ),
            parent_spatial_structure_draws=(
                parent_structure_draws
            ),
            campaign_structure_weights=(
                campaign_structure_weights
            ),
        )

        observed_conditional = family_result[
            "observed_conditional_proportions"
        ]

        summary_equal = posterior_summary(
            draws=equal_draws,
            mode_labels=family_result[
                "subfamily_labels"
            ],
            observed_proportions=(
                observed_conditional
            ),
            support_threshold=(
                support_threshold
            ),
            analysis_name=(
                f"{family_label}_conditional_equal"
            ),
        )

        summary_adjusted = posterior_summary(
            draws=adjusted_draws,
            mode_labels=family_result[
                "subfamily_labels"
            ],
            observed_proportions=(
                observed_conditional
            ),
            support_threshold=(
                support_threshold
            ),
            analysis_name=(
                f"{family_label}_conditional_adjusted"
            ),
        )

        output[
            family_label
        ] = {
            **family_result,
            "equal_draws": equal_draws,
            "adjusted_draws": adjusted_draws,
            "conditional_weights": (
                conditional_weights
            ),
            "summary_equal": (
                summary_equal
            ),
            "summary_adjusted": (
                summary_adjusted
            ),
        }

    return output

def save_conditional_support_plot(
    conditional_results: dict,
    output_file: Path,
):
    """
    Plot P(subfamily | spatial family) posterior distributions.

    One panel is generated per parent spatial family.
    """
    family_labels = sorted(
        conditional_results,
        key=lambda label: int(
            label[
                1:
            ]
        ),
    )

    n_families = len(
        family_labels
    )

    n_columns = 2

    n_rows = int(
        math.ceil(
            n_families
            / n_columns
        )
    )

    figure, axes = plt.subplots(
        n_rows,
        n_columns,
        figsize=(
            14.0,
            max(
                4.5,
                4.2
                * n_rows,
            ),
        ),
        squeeze=False,
    )

    flat_axes = axes.ravel()

    for panel_index, family_label in enumerate(
        family_labels
    ):
        axis = flat_axes[
            panel_index
        ]

        result = conditional_results[
            family_label
        ]

        labels = result[
            "subfamily_labels"
        ]

        equal_draws = result[
            "equal_draws"
        ]

        adjusted_draws = result[
            "adjusted_draws"
        ]

        adjusted_summary = (
            result[
                "summary_adjusted"
            ]
            .set_index(
                "mode"
            )
        )

        ordered_labels = (
            result[
                "summary_adjusted"
            ][
                "mode"
            ]
            .tolist()
        )

        label_to_index = {
            label: index
            for index, label in enumerate(
                labels
            )
        }

        positions = np.arange(
            len(
                ordered_labels
            )
        ) * 2.0

        equal_positions = positions
        adjusted_positions = (
            positions
            + 0.65
        )

        equal_values = [
            equal_draws[
                :,
                label_to_index[
                    label
                ],
            ]
            for label in ordered_labels
        ]

        adjusted_values = [
            adjusted_draws[
                :,
                label_to_index[
                    label
                ],
            ]
            for label in ordered_labels
        ]

        equal_violin = axis.violinplot(
            equal_values,
            positions=equal_positions,
            vert=False,
            widths=0.55,
            showmeans=False,
            showmedians=True,
            showextrema=False,
        )

        adjusted_violin = axis.violinplot(
            adjusted_values,
            positions=adjusted_positions,
            vert=False,
            widths=0.55,
            showmeans=False,
            showmedians=True,
            showextrema=False,
        )

        for body in equal_violin[
            "bodies"
        ]:
            body.set_facecolor(
                "#7aa6c2"
            )
            body.set_alpha(
                0.70
            )

        for body in adjusted_violin[
            "bodies"
        ]:
            body.set_facecolor(
                "#d98c4c"
            )
            body.set_alpha(
                0.70
            )

        axis.set_yticks(
            positions
            + 0.325
        )

        axis.set_yticklabels(
            ordered_labels
        )

        axis.set_xlim(
            0.0,
            1.0,
        )

        axis.set_xlabel(
            (
                "Conditional posterior support "
                f"P(subfamily | {family_label})"
            )
        )

        axis.set_title(
            (
                f"{family_label}: orientation/contact "
                "resolution conditional on position"
            )
        )

        # Annotate adjusted posterior medians.
        for position, label in zip(
            positions,
            ordered_labels,
        ):
            median = float(
                adjusted_summary.loc[
                    label,
                    "posterior_median",
                ]
            )

            axis.text(
                min(
                    median + 0.02,
                    0.95,
                ),
                position + 0.325,
                f"{median:.2f}",
                va="center",
                fontsize=8,
            )

    for unused_index in range(
        n_families,
        len(
            flat_axes
        ),
    ):
        flat_axes[
            unused_index
        ].axis(
            "off"
        )

    figure.suptitle(
        (
            "Conditional posterior support for "
            "orientation/contact modes"
        ),
        fontsize=14,
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

def build_hierarchical_support_summary(
    spatial_adjusted_summary: pd.DataFrame,
    conditional_results: dict,
    joint_adjusted_summary: pd.DataFrame,
) -> pd.DataFrame:
    spatial_lookup = (
        spatial_adjusted_summary
        .set_index(
            "mode"
        )
    )

    joint_lookup = (
        joint_adjusted_summary
        .set_index(
            "mode"
        )
    )

    rows = []

    for family_label, result in (
        conditional_results.items()
    ):
        conditional_summary = (
            result[
                "summary_adjusted"
            ]
        )

        for row in conditional_summary.itertuples():
            subfamily_label = row.mode

            rows.append(
                {
                    "spatial_family": (
                        family_label
                    ),
                    "subfamily": (
                        subfamily_label
                    ),

                    "spatial_support_median": float(
                        spatial_lookup.loc[
                            family_label,
                            "posterior_median",
                        ]
                    ),
                    "spatial_support_ci95_lower": float(
                        spatial_lookup.loc[
                            family_label,
                            "ci95_lower",
                        ]
                    ),
                    "spatial_support_ci95_upper": float(
                        spatial_lookup.loc[
                            family_label,
                            "ci95_upper",
                        ]
                    ),

                    "conditional_support_median": float(
                        row.posterior_median
                    ),
                    "conditional_support_ci95_lower": float(
                        row.ci95_lower
                    ),
                    "conditional_support_ci95_upper": float(
                        row.ci95_upper
                    ),

                    "joint_support_median": float(
                        joint_lookup.loc[
                            subfamily_label,
                            "posterior_median",
                        ]
                    ),
                    "joint_support_ci95_lower": float(
                        joint_lookup.loc[
                            subfamily_label,
                            "ci95_lower",
                        ]
                    ),
                    "joint_support_ci95_upper": float(
                        joint_lookup.loc[
                            subfamily_label,
                            "ci95_upper",
                        ]
                    ),

                    "conditional_probability_rank_1": float(
                        row.probability_rank_1
                    ),
                }
            )

    return (
        pd.DataFrame(
            rows
        )
        .sort_values(
            [
                "spatial_support_median",
                "conditional_support_median",
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