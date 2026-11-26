from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt


# ============================================================
# Fingerprint parsing
# ============================================================

def parse_contact_fingerprint(
    fingerprint_text: str,
) -> set[str]:
    """
    Parse the semicolon-delimited fingerprint stored in pose_table.

    Expected features include:

        R0123:contact
        R0123:polar_candidate
    """
    if pd.isna(
        fingerprint_text
    ):
        return set()

    return {
        feature.strip()
        for feature in str(
            fingerprint_text
        ).split(
            ";"
        )
        if feature.strip()
    }


def split_contact_feature(
    feature: str,
) -> tuple[str, str]:
    """
    Split:

        R0123:polar_candidate

    into:

        ("R0123", "polar_candidate")
    """
    if ":" not in feature:
        return (
            feature,
            "contact",
        )

    residue_label, contact_type = (
        feature.split(
            ":",
            1,
        )
    )

    return (
        residue_label,
        contact_type,
    )


def sequence_position_from_label(
    residue_label: str,
) -> int:
    """
    Convert R0123 to one-based reference-sequence position 123.
    """
    match = re.fullmatch(
        r"R(\d+)",
        residue_label,
    )

    if not match:
        raise ValueError(
            f"Invalid residue label: {residue_label}"
        )

    return int(
        match.group(
            1
        )
    )


# ============================================================
# Pose-level contact table
# ============================================================

def build_pose_contact_table(
    pose_table: pd.DataFrame,
) -> pd.DataFrame:
    required_columns = {
        "source_name",
        "pose_number",
        "global_spatial_family_id",
        "orientation_contact_subfamily_id",
        "contact_fingerprint",
    }

    missing = (
        required_columns
        - set(
            pose_table.columns
        )
    )

    if missing:
        raise ValueError(
            "Pocket-evidence analysis is missing pose-table "
            f"columns: {sorted(missing)}"
        )

    rows = []

    for pose_row in pose_table.itertuples(
        index=False
    ):
        mode_label = (
            f"F{int(pose_row.global_spatial_family_id)}."
            f"{int(pose_row.orientation_contact_subfamily_id)}"
        )

        features = parse_contact_fingerprint(
            pose_row.contact_fingerprint
        )

        for feature in features:
            (
                residue_label,
                contact_type,
            ) = split_contact_feature(
                feature
            )

            rows.append(
                {
                    "source_name": (
                        pose_row.source_name
                    ),
                    "pose_number": int(
                        pose_row.pose_number
                    ),
                    "spatial_family": (
                        f"F{int(pose_row.global_spatial_family_id)}"
                    ),
                    "mode": mode_label,
                    "residue_label": residue_label,
                    "reference_sequence_position": (
                        sequence_position_from_label(
                            residue_label
                        )
                    ),
                    "contact_type": contact_type,
                }
            )

    return pd.DataFrame(
        rows
    )


# ============================================================
# Mode-level contact probabilities
# ============================================================

def calculate_mode_contact_probabilities(
    pose_table: pd.DataFrame,
    pose_contacts: pd.DataFrame,
) -> pd.DataFrame:
    """
    Calculate P(contact with residue/type | mode) from poses assigned
    to each orientation/contact mode.
    """
    pose_modes = pose_table[
        [
            "source_name",
            "pose_number",
            "global_spatial_family_id",
            "orientation_contact_subfamily_id",
        ]
    ].copy()

    pose_modes[
        "mode"
    ] = [
        (
            f"F{int(family_id)}."
            f"{int(subfamily_id)}"
        )
        for family_id, subfamily_id in zip(
            pose_modes[
                "global_spatial_family_id"
            ],
            pose_modes[
                "orientation_contact_subfamily_id"
            ],
        )
    ]

    mode_pose_counts = (
        pose_modes
        .groupby(
            "mode"
        )
        .size()
        .rename(
            "n_mode_poses"
        )
    )

    if pose_contacts.empty:
        return pd.DataFrame(
            columns=[
                "mode",
                "residue_label",
                "reference_sequence_position",
                "contact_type",
                "n_mode_poses",
                "n_contacting_poses",
                "contact_probability_given_mode",
            ]
        )

    contact_counts = (
        pose_contacts
        .drop_duplicates(
            [
                "source_name",
                "pose_number",
                "mode",
                "residue_label",
                "contact_type",
            ]
        )
        .groupby(
            [
                "mode",
                "residue_label",
                "reference_sequence_position",
                "contact_type",
            ]
        )
        .size()
        .rename(
            "n_contacting_poses"
        )
        .reset_index()
    )

    contact_counts[
        "n_mode_poses"
    ] = contact_counts[
        "mode"
    ].map(
        mode_pose_counts
    )

    contact_counts[
        "contact_probability_given_mode"
    ] = (
        contact_counts[
            "n_contacting_poses"
        ]
        / contact_counts[
            "n_mode_poses"
        ]
    )

    return contact_counts


# ============================================================
# Structure-level contact probabilities
# ============================================================

def calculate_structure_contact_probabilities(
    pose_table: pd.DataFrame,
    pose_contacts: pd.DataFrame,
) -> pd.DataFrame:
    """
    Calculate residue-contact probability separately within each
    structure-level ensemble.
    """
    structure_pose_counts = (
        pose_table
        .groupby(
            "source_name"
        )
        .size()
        .rename(
            "n_structure_poses"
        )
    )

    contact_rows = (
        pose_contacts[
            pose_contacts[
                "contact_type"
            ] == "contact"
        ]
        .drop_duplicates(
            [
                "source_name",
                "pose_number",
                "residue_label",
            ]
        )
    )

    counts = (
        contact_rows
        .groupby(
            [
                "source_name",
                "residue_label",
                "reference_sequence_position",
            ]
        )
        .size()
        .rename(
            "n_contacting_poses"
        )
        .reset_index()
    )

    counts[
        "n_structure_poses"
    ] = counts[
        "source_name"
    ].map(
        structure_pose_counts
    )

    counts[
        "structure_contact_probability"
    ] = (
        counts[
            "n_contacting_poses"
        ]
        / counts[
            "n_structure_poses"
        ]
    )

    return counts


# ============================================================
# Posterior-weighted residue evidence
# ============================================================

def _adjusted_joint_mode_summary(
    joint_summary: pd.DataFrame,
) -> pd.DataFrame:
    adjusted = joint_summary[
        joint_summary[
            "analysis"
        ] == "joint_similarity_adjusted"
    ].copy()

    if adjusted.empty:
        raise RuntimeError(
            "No joint_similarity_adjusted rows were found in "
            "the subfamily posterior summary."
        )

    return adjusted


def calculate_posterior_weighted_residue_evidence(
    mode_contacts: pd.DataFrame,
    joint_summary: pd.DataFrame,
    structure_contacts: pd.DataFrame,
    total_structures: int,
) -> pd.DataFrame:
    """
    Combine:

        campaign P(mode)
        x
        P(residue contact | mode)

    yielding a posterior-weighted residue contact score.

    This score approximates the share of campaign-supported joint-mode
    probability associated with contacting a residue.
    """
    adjusted_modes = _adjusted_joint_mode_summary(
        joint_summary
    )

    mode_support = adjusted_modes.set_index(
        "mode"
    )[
        [
            "posterior_mean",
            "posterior_median",
            "ci95_lower",
            "ci95_upper",
        ]
    ]

    any_contacts = mode_contacts[
        mode_contacts[
            "contact_type"
        ] == "contact"
    ].copy()

    polar_contacts = mode_contacts[
        mode_contacts[
            "contact_type"
        ] == "polar_candidate"
    ].copy()

    any_contacts[
        "mode_support_mean"
    ] = any_contacts[
        "mode"
    ].map(
        mode_support[
            "posterior_mean"
        ]
    )

    any_contacts[
        "mode_support_median"
    ] = any_contacts[
        "mode"
    ].map(
        mode_support[
            "posterior_median"
        ]
    )

    any_contacts[
        "weighted_contact_mean"
    ] = (
        any_contacts[
            "mode_support_mean"
        ]
        * any_contacts[
            "contact_probability_given_mode"
        ]
    )

    any_contacts[
        "weighted_contact_median"
    ] = (
        any_contacts[
            "mode_support_median"
        ]
        * any_contacts[
            "contact_probability_given_mode"
        ]
    )

    residue_summary = (
        any_contacts
        .groupby(
            [
                "residue_label",
                "reference_sequence_position",
            ]
        )
        .agg(
            posterior_weighted_contact_mean=(
                "weighted_contact_mean",
                "sum",
            ),
            posterior_weighted_contact_median=(
                "weighted_contact_median",
                "sum",
            ),
            maximum_mode_contact_probability=(
                "contact_probability_given_mode",
                "max",
            ),
            n_contacting_modes=(
                "mode",
                "nunique",
            ),
        )
        .reset_index()
    )

    # Posterior-weighted polar-candidate evidence.
    if not polar_contacts.empty:
        polar_contacts[
            "mode_support_mean"
        ] = polar_contacts[
            "mode"
        ].map(
            mode_support[
                "posterior_mean"
            ]
        )

        polar_contacts[
            "weighted_polar_mean"
        ] = (
            polar_contacts[
                "mode_support_mean"
            ]
            * polar_contacts[
                "contact_probability_given_mode"
            ]
        )

        polar_summary = (
            polar_contacts
            .groupby(
                [
                    "residue_label",
                    "reference_sequence_position",
                ]
            )
            .agg(
                posterior_weighted_polar_candidate=(
                    "weighted_polar_mean",
                    "sum",
                )
            )
            .reset_index()
        )

        residue_summary = residue_summary.merge(
            polar_summary,
            on=[
                "residue_label",
                "reference_sequence_position",
            ],
            how="left",
        )

    else:
        residue_summary[
            "posterior_weighted_polar_candidate"
        ] = 0.0

    residue_summary[
        "posterior_weighted_polar_candidate"
    ] = residue_summary[
        "posterior_weighted_polar_candidate"
    ].fillna(
        0.0
    )

    # Structure-level replication and variability.
    structure_pivot = (
        structure_contacts
        .pivot(
            index="residue_label",
            columns="source_name",
            values="structure_contact_probability",
        )
        .fillna(
            0.0
        )
    )

    structure_metrics = []

    for residue_label, values in (
        structure_pivot.iterrows()
    ):
        probabilities = values.to_numpy(
            dtype=float
        )

        mean_value = float(
            np.mean(
                probabilities
            )
        )

        standard_deviation = float(
            np.std(
                probabilities
            )
        )

        coefficient_of_variation = (
            standard_deviation
            / mean_value
            if mean_value > 0.0
            else np.nan
        )

        structure_metrics.append(
            {
                "residue_label": residue_label,
                "mean_structure_contact_probability": (
                    mean_value
                ),
                "structure_contact_sd": (
                    standard_deviation
                ),
                "structure_contact_cv": (
                    coefficient_of_variation
                ),
                "n_supporting_structures_10pct": int(
                    np.sum(
                        probabilities >= 0.10
                    )
                ),
                "structure_support_fraction_10pct": float(
                    np.sum(
                        probabilities >= 0.10
                    )
                    / total_structures
                ),
            }
        )

    structure_metrics_df = pd.DataFrame(
        structure_metrics
    )

    residue_summary = residue_summary.merge(
        structure_metrics_df,
        on="residue_label",
        how="left",
    )

    # Mode specificity:
    # maximum fraction of a residue's weighted evidence coming from
    # any one mode.
    mode_contributions = (
        any_contacts
        .groupby(
            [
                "residue_label",
                "mode",
            ]
        )[
            "weighted_contact_mean"
        ]
        .sum()
        .reset_index()
    )

    specificity_rows = []

    for residue_label, group in (
        mode_contributions.groupby(
            "residue_label"
        )
    ):
        contributions = group[
            "weighted_contact_mean"
        ].to_numpy(
            dtype=float
        )

        total = float(
            contributions.sum()
        )

        if total > 0.0:
            mode_specificity = float(
                contributions.max()
                / total
            )

            dominant_mode = group.iloc[
                int(
                    np.argmax(
                        contributions
                    )
                )
            ][
                "mode"
            ]

        else:
            mode_specificity = np.nan
            dominant_mode = ""

        specificity_rows.append(
            {
                "residue_label": residue_label,
                "mode_specificity": mode_specificity,
                "dominant_contact_mode": dominant_mode,
            }
        )

    residue_summary = residue_summary.merge(
        pd.DataFrame(
            specificity_rows
        ),
        on="residue_label",
        how="left",
    )

    return residue_summary


# ============================================================
# Residue classification
# ============================================================

def classify_residue_evidence(
    residue_table: pd.DataFrame,
    conserved_support_threshold: float,
    moderate_support_threshold: float,
    maximum_conserved_cv: float,
    minimum_supporting_structures: int,
    mode_specificity_threshold: float,
    weak_support_threshold: float,
) -> pd.DataFrame:
    """
    Classify the strength and reproducibility of residue-engagement
    evidence.

    These categories describe structural evidence in the current
    pose ensembles. They do not predict affinity contributions,
    energetic favourability, or whether targeting a residue will
    improve a compound.
    """
    result = residue_table.copy()

    categories = []
    rationales = []

    for row in result.itertuples(
        index=False
    ):
        weighted_support = float(
            row.posterior_weighted_contact_mean
        )

        supporting_structures = int(
            row.n_supporting_structures_10pct
        )

        cv = (
            float(
                row.structure_contact_cv
            )
            if np.isfinite(
                row.structure_contact_cv
            )
            else np.inf
        )

        specificity = (
            float(
                row.mode_specificity
            )
            if np.isfinite(
                row.mode_specificity
            )
            else 0.0
        )

        if (
            weighted_support
            >= conserved_support_threshold
            and supporting_structures
            >= minimum_supporting_structures
            and cv
            <= maximum_conserved_cv
        ):
            category = (
                "conserved_engagement"
            )

            rationale = (
                "High posterior-weighted contact support, "
                "replicated across structure-level ensembles, "
                "with relatively low between-structure "
                "variability."
            )

        elif (
            weighted_support
            >= moderate_support_threshold
            and specificity
            >= mode_specificity_threshold
        ):
            category = (
                "mode_specific_engagement"
            )

            rationale = (
                "Meaningful contact evidence is concentrated "
                "in one orientation/contact mode; interpretation "
                "therefore depends strongly on support for that "
                "mode."
            )

        elif (
            weighted_support
            >= moderate_support_threshold
            and supporting_structures
            >= minimum_supporting_structures
        ):
            category = (
                "variable_engagement"
            )

            rationale = (
                "Meaningful contact evidence is reproduced "
                "across structures, but its prevalence varies "
                "between structure-level ensembles or supported "
                "binding modes."
            )

        elif (
            weighted_support
            <= weak_support_threshold
        ):
            category = (
                "weakly_evidenced"
            )

            rationale = (
                "Little of the posterior-supported binding "
                "landscape contacts this residue. This does not "
                "establish that the residue is inaccessible or "
                "energetically unfavourable."
            )

        else:
            category = (
                "limited_or_context_dependent"
            )

            rationale = (
                "Contact evidence is present but limited, "
                "variable, insufficiently replicated, or below "
                "the criteria for the stronger evidence classes."
            )

        categories.append(
            category
        )

        rationales.append(
            rationale
        )

    result[
        "evidence_category"
    ] = categories

    result[
        "classification_rationale"
    ] = rationales

    result[
        "engagement_support"
    ] = result[
        "posterior_weighted_contact_mean"
    ]

    result[
        "between_structure_consistency"
    ] = (
        1.0
        / (
            1.0
            + result[
                "structure_contact_cv"
            ]
        )
    )

    return result


# ============================================================
# Residue-mode matrix
# ============================================================

def build_residue_mode_matrix(
    mode_contacts: pd.DataFrame,
    joint_summary: pd.DataFrame,
) -> pd.DataFrame:
    """
    Matrix entries are:

        posterior P(mode)
        x
        P(residue contact | mode)
    """
    adjusted_modes = _adjusted_joint_mode_summary(
        joint_summary
    )

    support_map = adjusted_modes.set_index(
        "mode"
    )[
        "posterior_mean"
    ]

    contacts = mode_contacts[
        mode_contacts[
            "contact_type"
        ] == "contact"
    ].copy()

    contacts[
        "weighted_contact"
    ] = (
        contacts[
            "contact_probability_given_mode"
        ]
        * contacts[
            "mode"
        ].map(
            support_map
        )
    )

    matrix = (
        contacts
        .pivot_table(
            index="residue_label",
            columns="mode",
            values="weighted_contact",
            aggfunc="sum",
            fill_value=0.0,
        )
    )

    matrix[
        "_total"
    ] = matrix.sum(
        axis=1
    )

    matrix = (
        matrix
        .sort_values(
            "_total",
            ascending=False,
        )
        .drop(
            columns=[
                "_total",
            ]
        )
    )

    column_support = matrix.sum(
        axis=0
    )

    matrix = matrix[
        column_support.sort_values(
            ascending=False
        ).index
    ]

    return matrix


# ============================================================
# Plots
# ============================================================

CATEGORY_ORDER = [
    "conserved_engagement",
    "variable_engagement",
    "mode_specific_engagement",
    "limited_or_context_dependent",
    "weakly_evidenced",
]

CATEGORY_LABELS = {
    "conserved_engagement": "Conserved engagement",
    "variable_engagement": "Variable engagement",
    "mode_specific_engagement": "Mode-specific engagement",
    "limited_or_context_dependent": "Limited / context-dependent",
    "weakly_evidenced": "Weakly evidenced",
}

CATEGORY_COLOURS = {
    "conserved_engagement": "#00A651",          # vivid green
    "variable_engagement": "#0077FF",           # vivid blue
    "mode_specific_engagement": "#FFB000",      # vivid amber
    "limited_or_context_dependent": "#8A2BE2",  # vivid purple
    "weakly_evidenced": "#A6A6A6",              # neutral grey
}

def save_residue_evidence_plot(
    residue_table: pd.DataFrame,
    output_file: Path,
    top_n: int,
):
    display = (
        residue_table
        .sort_values(
            "posterior_weighted_contact_mean",
            ascending=False,
        )
        .head(
            top_n
        )
        .sort_values(
            "posterior_weighted_contact_mean",
            ascending=True,
        )
    )

    colours = [
        CATEGORY_COLOURS.get(
            category,
            "#bdbdbd",
        )
        for category in display[
            "evidence_category"
        ]
    ]

    figure_height = max(
        6.0,
        0.35
        * len(
            display
        )
        + 2.0,
    )

    figure, axis = plt.subplots(
        figsize=(
            11.0,
            figure_height,
        )
    )

    positions = np.arange(
        len(
            display
        )
    )

    axis.barh(
        positions,
        display[
            "posterior_weighted_contact_mean"
        ],
        color=colours,
        alpha=0.85,
    )

    axis.set_yticks(
        positions
    )

    axis.set_yticklabels(
        display[
            "residue_label"
        ]
    )

    axis.set_xlabel(
        "Posterior-weighted residue contact support"
    )

    axis.set_ylabel(
        "Reference-sequence position"
    )

    axis.set_title(
        "Campaign structural evidence by pocket residue"
    )

    displayed_categories = set(
        display[
            "evidence_category"
        ]
    )

    legend_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="s",
            linestyle="",
            markerfacecolor=(
                CATEGORY_COLOURS[
                    category
                ]
            ),
            markeredgecolor="black",
            markeredgewidth=0.5,
            markersize=9,
            label=(
                CATEGORY_LABELS[
                    category
                ]
            ),
        )
        for category in CATEGORY_ORDER
        if category in displayed_categories
    ]

    axis.legend(
        handles=legend_handles,
        loc="best",
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


def save_residue_engagement_evidence_map(
    residue_table: pd.DataFrame,
    output_file: Path,
    conserved_support_threshold: float,
    maximum_conserved_cv: float,
):
    """
    Plot residue-level engagement evidence.

    x:
        Posterior-weighted residue contact support across the
        supported campaign binding-mode landscape.

    y:
        Consistency of contact prevalence across structure-level
        ensembles, calculated as 1 / (1 + CV).

    colour:
        Structural-evidence category. Categories describe support,
        replication, variability and mode dependence. They do not
        represent affinity contributions or energetic favourability.

    size:
        Posterior-weighted geometric polar-proximity evidence.
        This is not a validated hydrogen-bond probability.

    The dashed lines show only the support and consistency criteria
    that contribute to the conserved-engagement classification.
    Other categories also depend on replication and mode specificity.
    """
    table = residue_table.copy()

    consistency = (
        1.0
        / (
            1.0
            + table[
                "structure_contact_cv"
            ]
            .replace(
                [
                    np.inf,
                    -np.inf,
                ],
                np.nan,
            )
            .fillna(
                10.0
            )
        )
    )

    point_sizes = (
        40.0
        + 700.0
        * table[
            "posterior_weighted_polar_candidate"
        ].clip(
            lower=0.0
        )
    )

    conserved_consistency_threshold = (
        1.0
        / (
            1.0
            + maximum_conserved_cv
        )
    )

    figure, axis = plt.subplots(
        figsize=(
            10.0,
            7.5,
        )
    )

    # Draw categories in explicit evidence order so the legend
    # follows the same order.
    for category in CATEGORY_ORDER:
        group = table[
            table[
                "evidence_category"
            ] == category
        ]

        if group.empty:
            continue

        indices = group.index

        axis.scatter(
            group[
                "posterior_weighted_contact_mean"
            ],
            consistency.loc[
                indices
            ],
            s=point_sizes.loc[
                indices
            ],
            color=(
                CATEGORY_COLOURS[
                    category
                ]
            ),
            alpha=0.88,
            label=(
                CATEGORY_LABELS[
                    category
                ]
            ),
            edgecolor="black",
            linewidth=0.55,
            zorder=3,
        )

    # Label the highest-support residues.
    top_labels = (
        table
        .sort_values(
            "posterior_weighted_contact_mean",
            ascending=False,
        )
        .head(
            min(
                20,
                len(
                    table
                ),
            )
        )
    )

    for row in top_labels.itertuples():
        axis.annotate(
            row.residue_label,
            (
                row.posterior_weighted_contact_mean,
                consistency.loc[
                    row.Index
                ],
            ),
            xytext=(
                5,
                5,
            ),
            textcoords="offset points",
            fontsize=8,
            zorder=4,
        )

    # Conserved-engagement boundaries.
    axis.axvline(
        conserved_support_threshold,
        color="#4D4D4D",
        linestyle="--",
        linewidth=1.1,
        alpha=0.75,
        zorder=1,
    )

    axis.axhline(
        conserved_consistency_threshold,
        color="#4D4D4D",
        linestyle="--",
        linewidth=1.1,
        alpha=0.75,
        zorder=1,
    )

    # Label the vertical criterion near the bottom of the plot.
    axis.annotate(
        (
            "High support criterion\n"
            f"support >= {conserved_support_threshold:.2f}"
        ),
        xy=(
            conserved_support_threshold,
            0.02,
        ),
        xycoords=(
            "data",
            "axes fraction",
        ),
        xytext=(
            5,
            3,
        ),
        textcoords="offset points",
        ha="left",
        va="bottom",
        fontsize=8,
        color="#3D3D3D",
        rotation=90,
    )

    # Label the horizontal criterion on the right side.
    axis.annotate(
        (
            "High consistency criterion  "
            f"{conserved_consistency_threshold:.2f}"
        ),
        xy=(
            0.99,
            conserved_consistency_threshold,
        ),
        xycoords=(
            "axes fraction",
            "data",
        ),
        xytext=(
            -4,
            5,
        ),
        textcoords="offset points",
        ha="right",
        va="bottom",
        fontsize=8,
        color="#3D3D3D",
    )

    axis.set_xlabel(
        "Posterior-weighted engagement support"
    )

    axis.set_ylabel(
        "Cross-structure engagement consistency"
    )

    axis.set_ylim(
        0.0,
        1.05,
    )

    axis.set_title(
        "Residue engagement evidence across all cryo-em structures."
    )

    axis.legend(
        loc="best",
        title="Evidence classification",
        frameon=True,
    )

    # Method annotation. This is intentionally separate from the
    # threshold labels because the categories use more than x/y.
    figure.text(
        0.5,
        0.01,
        (
            "Consistency = 1 / (1 + between-structure CV). "
            "Dashed lines indicate high-support and high-consistency "
            "boundaries. Evidence classes additionally consider "
            "cross-structure replication and mode specificity."
        ),
        ha="center",
        va="bottom",
        fontsize=8,
        color="#444444",
    )

    figure.tight_layout(
        rect=(
            0.0,
            0.045,
            1.0,
            1.0,
        )
    )

    figure.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(
        figure
    )


def save_residue_mode_heatmap(
    matrix: pd.DataFrame,
    output_file: Path,
    top_n_residues: int,
    top_n_modes: int = 20,
):
    display = matrix.iloc[
        :top_n_residues,
        :top_n_modes,
    ]

    figure_width = max(
        9.0,
        0.55
        * len(
            display.columns
        )
        + 3.0,
    )

    figure_height = max(
        7.0,
        0.35
        * len(
            display.index
        )
        + 2.0,
    )

    figure, axis = plt.subplots(
        figsize=(
            figure_width,
            figure_height,
        )
    )

    image = axis.imshow(
        display.to_numpy(
            dtype=float
        ),
        aspect="auto",
        cmap="Blues",
        vmin=0.0,
        vmax=max(
            float(
                display.to_numpy(
                    dtype=float
                ).max()
            ),
            1.0e-6,
        ),
    )

    axis.set_xticks(
        np.arange(
            len(
                display.columns
            )
        )
    )

    axis.set_xticklabels(
        display.columns,
        rotation=90,
    )

    axis.set_yticks(
        np.arange(
            len(
                display.index
            )
        )
    )

    axis.set_yticklabels(
        display.index
    )

    axis.set_xlabel(
        "Orientation/contact mode"
    )

    axis.set_ylabel(
        "Reference-sequence position"
    )

    axis.set_title(
        "Posterior-weighted residue-contact evidence by mode"
    )

    colour_bar = figure.colorbar(
        image,
        ax=axis,
    )

    colour_bar.set_label(
        "P(mode) × P(contact | mode)"
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
# PyMOL residue mapping
# ============================================================

def write_reference_residue_mapping(
    alignment_reference_structure,
    output_file: Path,
    reference_chain_id: str,
) -> pd.DataFrame:
    """
    Map the canonical one-based sequence positions used in contact
    fingerprints onto residue identifiers in the aligned reference
    receptor.

    Only the sequence-selected reference receptor chain is used.
    This ensures:

        R0001, R0002, ...

    have exactly the same meaning in the fingerprints, CSV outputs
    and PyMOL selections.

    Parameters
    ----------
    alignment_reference_structure
        Receptor-only aligned reference Structure.

    output_file
        CSV output path.

    reference_chain_id
        Chain selected by sequence matching for the campaign reference
        receptor. The chain name is local to the reference structure
        and does not need to match chains in other structures.
    """
    output_file = Path(
        output_file
    )

    expected_columns = [
        "residue_label",
        "reference_sequence_position",
        "chain",
        "pdb_resnum",
        "insertion_code",
        "residue_name",
        "ca_atom_index",
    ]

    rows = []
    seen_residues = set()

    reference_chain_id = (
        reference_chain_id
        or ""
    ).strip()

    for atom in alignment_reference_structure.atom:
        atom_name = (
            getattr(
                atom,
                "pdbname",
                "",
            )
            or ""
        ).strip().upper()

        if atom_name != "CA":
            continue

        atom_chain = (
            getattr(
                atom,
                "chain",
                "",
            )
            or ""
        ).strip()

        if atom_chain != reference_chain_id:
            continue

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
            atom_chain,
            residue_number,
            insertion_code,
        )

        if residue_key in seen_residues:
            continue

        seen_residues.add(
            residue_key
        )

        sequence_position = len(
            rows
        ) + 1

        rows.append(
            {
                "residue_label": (
                    f"R{sequence_position:04d}"
                ),
                "reference_sequence_position": (
                    sequence_position
                ),
                "chain": atom_chain,
                "pdb_resnum": residue_number,
                "insertion_code": insertion_code,
                "residue_name": residue_name,
                "ca_atom_index": int(
                    atom.index
                ),
            }
        )

    # Supplying columns explicitly ensures an empty DataFrame still has
    # a defined schema, although an empty result is treated as an error.
    mapping = pd.DataFrame(
        rows,
        columns=expected_columns,
    )

    if mapping.empty:
        available_ca_chains = sorted(
            {
                (
                    getattr(
                        atom,
                        "chain",
                        "",
                    )
                    or ""
                ).strip()
                for atom in alignment_reference_structure.atom
                if (
                    (
                        getattr(
                            atom,
                            "pdbname",
                            "",
                        )
                        or ""
                    )
                    .strip()
                    .upper()
                    == "CA"
                )
            }
        )

        raise RuntimeError(
            "No C-alpha residues were found for the selected "
            f"reference receptor chain {reference_chain_id!r}. "
            f"Available C-alpha chains: {available_ca_chains}"
        )

    if mapping[
        "residue_label"
    ].duplicated().any():
        raise RuntimeError(
            "Reference residue mapping contains duplicate "
            "canonical residue labels."
        )

    if mapping[
        [
            "chain",
            "pdb_resnum",
            "insertion_code",
        ]
    ].duplicated().any():
        duplicates = mapping[
            mapping[
                [
                    "chain",
                    "pdb_resnum",
                    "insertion_code",
                ]
            ].duplicated(
                keep=False,
            )
        ]

        raise RuntimeError(
            "Reference residue mapping contains duplicate PDB "
            "residue identities:\n"
            f"{duplicates.to_string(index=False)}"
        )

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    mapping.to_csv(
        output_file,
        index=False,
    )

    return mapping


def write_pocket_evidence_pml(
    residue_table: pd.DataFrame,
    residue_mapping: pd.DataFrame,
    receptor_pdb_name: str,
    output_file: Path,
):
    """
    Generate a portable PML using a receptor PDB beside the script.
    """
    required_residue_columns = {
        "residue_label",
        "reference_sequence_position",
        "evidence_category",
    }

    required_mapping_columns = {
        "residue_label",
        "reference_sequence_position",
        "chain",
        "pdb_resnum",
        "insertion_code",
        "residue_name",
    }

    missing_residue_columns = (
        required_residue_columns
        - set(
            residue_table.columns
        )
    )

    missing_mapping_columns = (
        required_mapping_columns
        - set(
            residue_mapping.columns
        )
    )

    if missing_residue_columns:
        raise RuntimeError(
            "Residue evidence table is missing columns required "
            f"for PyMOL export: {sorted(missing_residue_columns)}"
        )

    if missing_mapping_columns:
        raise RuntimeError(
            "Reference residue mapping is missing columns required "
            f"for PyMOL export: {sorted(missing_mapping_columns)}. "
            f"Available columns: {list(residue_mapping.columns)}"
        )

    merged = residue_table.merge(
        residue_mapping,
        on=[
            "residue_label",
            "reference_sequence_position",
        ],
        how="inner",
        validate="one_to_one",
    )

    if merged.empty:
        evidence_labels = sorted(
            residue_table[
                "residue_label"
            ].astype(
                str
            ).unique()
            .tolist()
        )

        mapping_labels = sorted(
            residue_mapping[
                "residue_label"
            ].astype(
                str
            ).unique()
            .tolist()
        )

        raise RuntimeError(
            "Pocket evidence and reference residue mapping have "
            "no matching residue labels. "
            f"First evidence labels: {evidence_labels[:10]}; "
            f"first mapping labels: {mapping_labels[:10]}"
        )

    lines = [
        "reinitialize",
        "bg_color white",
        (
            f"load {receptor_pdb_name}, "
            "receptor"
        ),
        "hide everything, receptor",
        "show cartoon, receptor",
        "color gray80, receptor",
    ]

    for category, colour in (
        CATEGORY_COLOURS.items()
    ):
        category_rows = merged[
            merged[
                "evidence_category"
            ] == category
        ]

        selections = []

        for row in category_rows.itertuples(
            index=False
        ):
            clause = (
                f"(chain {row.chain} "
                f"and resi {row.pdb_resnum}"
            )

            if row.insertion_code:
                clause += (
                    f"{row.insertion_code}"
                )

            clause += ")"

            selections.append(
                clause
            )

        if not selections:
            continue

        selection_name = (
            f"pocket_{category}"
        )

        selection_expression = (
            " or ".join(
                selections
            )
        )

        lines.extend(
            [
                (
                    f"select {selection_name}, "
                    f"({selection_expression})"
                ),
                (
                    f"show sticks, "
                    f"{selection_name}"
                ),
                (
                    f"color {colour}, "
                    f"({selection_name} and elem C)"
                ),
                (
                    f"color blue, "
                    f"({selection_name} and elem N)"
                ),
                (
                    f"color red, "
                    f"({selection_name} and elem O)"
                ),
                (
                    f"color yellow, "
                    f"({selection_name} and elem S)"
                ),
            ]
        )

    lines.extend(
        [
            "set stick_radius, 0.18",
            "orient receptor",
            "zoom receptor, 8",
        ]
    )

    output_file.write_text(
        "\n".join(
            lines
        )
        + "\n"
    )


# ============================================================
# Main entry point
# ============================================================

def run_pocket_evidence_analysis(
    pose_table: pd.DataFrame,
    joint_summary: pd.DataFrame,
    reference_receptor,
    reference_chain_id: str,
    reference_receptor_pdb_name: str,
    output_dir: Path,
    config: dict,
) -> dict:
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

    conserved_support_threshold = float(
        config.get(
            "conserved_support_threshold",
            0.35,
        )
    )

    moderate_support_threshold = float(
        config.get(
            "moderate_support_threshold",
            0.15,
        )
    )

    maximum_conserved_cv = float(
        config.get(
            "maximum_conserved_cv",
            0.50,
        )
    )

    minimum_supporting_structures = int(
        config.get(
            "minimum_supporting_structures",
            2,
        )
    )

    mode_specificity_threshold = float(
        config.get(
            "mode_specificity_threshold",
            0.60,
        )
    )

    weak_support_threshold = float(
        config.get(
            "weak_support_threshold",
            0.05,
        )
    )

    top_n_residues_plot = int(
        config.get(
            "top_n_residues_plot",
            40,
        )
    )

    pose_contacts = build_pose_contact_table(
        pose_table
    )

    pose_contacts.to_csv(
        output_dir
        / "pose_contact_features.csv",
        index=False,
    )

    mode_contacts = calculate_mode_contact_probabilities(
        pose_table=pose_table,
        pose_contacts=pose_contacts,
    )

    mode_contacts.to_csv(
        output_dir
        / "mode_contact_probabilities.csv",
        index=False,
    )

    structure_contacts = (
        calculate_structure_contact_probabilities(
            pose_table=pose_table,
            pose_contacts=pose_contacts,
        )
    )

    structure_contacts.to_csv(
        output_dir
        / "structure_contact_probabilities.csv",
        index=False,
    )

    residue_evidence = (
        calculate_posterior_weighted_residue_evidence(
            mode_contacts=mode_contacts,
            joint_summary=joint_summary,
            structure_contacts=(
                structure_contacts
            ),
            total_structures=int(
                pose_table[
                    "source_name"
                ].nunique()
            ),
        )
    )

    residue_evidence = classify_residue_evidence(
        residue_table=residue_evidence,
        conserved_support_threshold=(
            conserved_support_threshold
        ),
        moderate_support_threshold=(
            moderate_support_threshold
        ),
        maximum_conserved_cv=(
            maximum_conserved_cv
        ),
        minimum_supporting_structures=(
            minimum_supporting_structures
        ),
        mode_specificity_threshold=(
            mode_specificity_threshold
        ),
        weak_support_threshold=(
            weak_support_threshold
        ),
    )

    residue_evidence = residue_evidence.sort_values(
        "posterior_weighted_contact_mean",
        ascending=False,
    ).reset_index(
        drop=True
    )

    residue_evidence.to_csv(
        output_dir
        / "pocket_residue_evidence.csv",
        index=False,
    )

    residue_mode_matrix = build_residue_mode_matrix(
        mode_contacts=mode_contacts,
        joint_summary=joint_summary,
    )

    residue_mode_matrix.to_csv(
        output_dir
        / "residue_mode_evidence_matrix.csv"
    )

    save_residue_evidence_plot(
        residue_table=residue_evidence,
        output_file=(
            output_dir
            / "residue_evidence_ranking.png"
        ),
        top_n=top_n_residues_plot,
    )

    save_residue_engagement_evidence_map(
        residue_table=residue_evidence,
        output_file=(
            output_dir
            / "residue_engagement_evidence_map.png"
        ),
        conserved_support_threshold=(
            conserved_support_threshold
        ),
        maximum_conserved_cv=(
            maximum_conserved_cv
        ),
    )

    save_residue_mode_heatmap(
        matrix=residue_mode_matrix,
        output_file=(
            output_dir
            / "residue_mode_heatmap.png"
        ),
        top_n_residues=top_n_residues_plot,
    )

    residue_mapping = write_reference_residue_mapping(
        alignment_reference_structure=(
            reference_receptor
        ),
        reference_chain_id=(
            reference_chain_id
        ),
        output_file=(
            output_dir
            / "reference_residue_mapping.csv"
        ),
    )

    pymol_script = (
        output_dir
        / "residue_engagement_evidence.pml"
    )

    if bool(
        config.get(
            "export_pymol",
            True,
        )
    ):
        write_pocket_evidence_pml(
            residue_table=residue_evidence,
            residue_mapping=residue_mapping,
            receptor_pdb_name=(
                reference_receptor_pdb_name
            ),
            output_file=pymol_script,
        )

    return {
        "residue_evidence": (
            output_dir
            / "pocket_residue_evidence.csv"
        ),
        "mode_contacts": (
            output_dir
            / "mode_contact_probabilities.csv"
        ),
        "engagement_evidence_plot": (
            output_dir
            / "residue_engagement_evidence_map.png"
        ),
        "residue_mode_heatmap": (
            output_dir
            / "residue_mode_heatmap.png"
        ),
        "pymol_script": (
            pymol_script
            if pymol_script.is_file()
            else None
        ),
    }