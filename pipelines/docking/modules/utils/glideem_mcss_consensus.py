from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from schrodinger import structure

from modules.utils.glideem_helpers import (
    compound_id,
    glide_score,
    smiles_value,
)

from modules.utils.mcss_docking_helpers import (
    export_consensus_compound_pose,
    export_combined_consensus_poses,
)

from pipeline.logger import setup_logger


logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True,
)


def collect_mcss_docking_scores(
    complex_jobs,
    poseviewers,
) -> pd.DataFrame:
    """
    Read every available MCSS PoseViewer and retain the best
    GlideScore for each compound within each individual MCSS dock.

    This function does not use experimental activity.
    """
    jobs_by_name = {
        job[
            "name"
        ]: job
        for job in complex_jobs
    }

    rows = []

    for dock_name, poseviewer in (
        poseviewers.items()
    ):
        if dock_name not in jobs_by_name:
            logger.warning(
                "[%s] No complex-job metadata; excluded from "
                "MCSS consensus analysis",
                dock_name,
            )
            continue

        job = jobs_by_name[
            dock_name
        ]

        poseviewer = (
            Path(
                poseviewer
            )
            .expanduser()
            .resolve()
        )

        if not poseviewer.is_file():
            logger.warning(
                "[%s] PoseViewer is missing; excluded from "
                "MCSS consensus analysis: %s",
                dock_name,
                poseviewer,
            )
            continue

        best_by_compound = {}
        pose_count_by_compound = {}

        for ct_index, st in enumerate(
            structure.StructureReader(
                str(
                    poseviewer
                )
            ),
            start=1,
        ):
            if st.property.get(
                "b_glide_receptor",
                False,
            ):
                continue

            cid = compound_id(
                st
            )

            score = glide_score(
                st
            )

            if (
                cid is None
                or score is None
            ):
                continue

            cid = str(
                cid
            )

            score = float(
                score
            )

            pose_count_by_compound[
                cid
            ] = (
                pose_count_by_compound.get(
                    cid,
                    0,
                )
                + 1
            )

            compound_pose_index = (
                pose_count_by_compound[
                    cid
                ]
            )

            existing = best_by_compound.get(
                cid
            )

            if (
                existing is None
                or score
                < existing[
                    "glide_score"
                ]
            ):
                best_by_compound[
                    cid
                ] = {
                    "dock_name": dock_name,
                    "source_name": job[
                        "source_name"
                    ],
                    "glideem_pose_number": int(
                        job[
                            "pose_number"
                        ]
                    ),
                    "compound_id": cid,
                    "smiles": smiles_value(
                        st
                    ),
                    "glide_score": score,
                    "compound_pose_index": int(
                        compound_pose_index
                    ),
                    "poseviewer_ct_index": int(
                        ct_index
                    ),
                    "poseviewer_path": str(
                        poseviewer
                    ),
                    "reference_ligand_path": str(
                        Path(
                            job[
                                "reference"
                            ]
                        )
                        .expanduser()
                        .resolve()
                    ),
                }

        rows.extend(
            best_by_compound.values()
        )

        logger.debug(
            "[%s] Collected best GlideScores for %d compounds",
            dock_name,
            len(
                best_by_compound
            ),
        )

    if not rows:
        raise RuntimeError(
            "No scored MCSS compounds were found in the "
            "available PoseViewers."
        )

    return pd.DataFrame(
        rows
    )


def calculate_mcss_consensus_ranking(
    dock_scores: pd.DataFrame,
    top_ks: tuple[int, ...] = (
        1,
        5,
        10,
    ),
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
]:
    """
    Rank compounds within each MCSS dock, aggregate across docks
    within each cryo-EM structure, and then aggregate structures
    with equal weight.
    """
    working = dock_scores.copy()

    working[
        "dock_rank"
    ] = (
        working
        .groupby(
            "dock_name"
        )[
            "glide_score"
        ]
        .rank(
            method="min",
            ascending=True,
        )
        .astype(
            int
        )
    )

    working[
        "n_compounds_in_dock"
    ] = (
        working
        .groupby(
            "dock_name"
        )[
            "compound_id"
        ]
        .transform(
            "nunique"
        )
        .astype(
            int
        )
    )

    denominator = (
        working[
            "n_compounds_in_dock"
        ]
        - 1
    )

    # 1.0 is the best compound in an individual MCSS dock.
    # 0.0 is the worst.
    working[
        "dock_percentile_score"
    ] = np.where(
        denominator > 0,
        1.0
        - (
            working[
                "dock_rank"
            ]
            - 1
        )
        / denominator,
        1.0,
    )

    working[
        "reciprocal_rank"
    ] = (
        1.0
        / working[
            "dock_rank"
        ]
    )

    for top_k in top_ks:
        working[
            f"is_top_{top_k}"
        ] = (
            working[
                "dock_rank"
            ]
            <= top_k
        )

    # --------------------------------------------------
    # First level: aggregate MCSS docks within one
    # cryo-EM source structure.
    # --------------------------------------------------

    structure_aggregations = {
        "structure_median_percentile": (
            "dock_percentile_score",
            "median",
        ),
        "structure_mean_percentile": (
            "dock_percentile_score",
            "mean",
        ),
        "structure_best_percentile": (
            "dock_percentile_score",
            "max",
        ),
        "structure_worst_percentile": (
            "dock_percentile_score",
            "min",
        ),
        "structure_mean_reciprocal_rank": (
            "reciprocal_rank",
            "mean",
        ),
        "structure_median_rank": (
            "dock_rank",
            "median",
        ),
        "structure_best_glide_score": (
            "glide_score",
            "min",
        ),
        "n_docks_evaluated": (
            "dock_name",
            "nunique",
        ),
    }

    for top_k in top_ks:
        structure_aggregations[
            f"structure_top_{top_k}_count"
        ] = (
            f"is_top_{top_k}",
            "sum",
        )

        structure_aggregations[
            f"structure_top_{top_k}_fraction"
        ] = (
            f"is_top_{top_k}",
            "mean",
        )

    structure_scores = (
        working
        .groupby(
            [
                "compound_id",
                "source_name",
            ],
            as_index=False,
        )
        .agg(
            **structure_aggregations
        )
    )

    # --------------------------------------------------
    # Second level: aggregate equally across cryo-EM
    # source structures.
    # --------------------------------------------------

    campaign_aggregations = {
        "consensus_score": (
            "structure_median_percentile",
            "mean",
        ),
        "median_structure_percentile": (
            "structure_median_percentile",
            "median",
        ),
        "minimum_structure_percentile": (
            "structure_median_percentile",
            "min",
        ),
        "maximum_structure_percentile": (
            "structure_median_percentile",
            "max",
        ),
        "mean_structure_reciprocal_rank": (
            "structure_mean_reciprocal_rank",
            "mean",
        ),
        "n_structures_evaluated": (
            "source_name",
            "nunique",
        ),
        "total_docks_evaluated": (
            "n_docks_evaluated",
            "sum",
        ),
    }

    for top_k in top_ks:
        campaign_aggregations[
            f"mean_structure_top_{top_k}_fraction"
        ] = (
            f"structure_top_{top_k}_fraction",
            "mean",
        )

        campaign_aggregations[
            f"total_top_{top_k}_count"
        ] = (
            f"structure_top_{top_k}_count",
            "sum",
        )

    consensus = (
        structure_scores
        .groupby(
            "compound_id",
            as_index=False,
        )
        .agg(
            **campaign_aggregations
        )
    )

    structure_iqr = (
        structure_scores
        .groupby(
            "compound_id"
        )[
            "structure_median_percentile"
        ]
        .quantile(
            0.75
        )
        - structure_scores
        .groupby(
            "compound_id"
        )[
            "structure_median_percentile"
        ]
        .quantile(
            0.25
        )
    )

    consensus[
        "structure_percentile_iqr"
    ] = consensus[
        "compound_id"
    ].map(
        structure_iqr
    ).fillna(
        0.0
    )

    consensus[
        "structure_rank_stability"
    ] = (
        1.0
        - consensus[
            "structure_percentile_iqr"
        ].clip(
            lower=0.0,
            upper=1.0,
        )
    )

    total_structures = int(
        working[
            "source_name"
        ].nunique()
    )

    total_docks = int(
        working[
            "dock_name"
        ].nunique()
    )

    consensus[
        "structure_coverage_fraction"
    ] = (
        consensus[
            "n_structures_evaluated"
        ]
        / total_structures
    )

    consensus[
        "dock_coverage_fraction"
    ] = (
        consensus[
            "total_docks_evaluated"
        ]
        / total_docks
    )

    # Add representative SMILES where available.
    smiles_lookup = (
        working
        .dropna(
            subset=[
                "smiles",
            ]
        )
        .drop_duplicates(
            subset=[
                "compound_id",
            ]
        )
        .set_index(
            "compound_id"
        )[
            "smiles"
        ]
    )

    consensus[
        "smiles"
    ] = consensus[
        "compound_id"
    ].map(
        smiles_lookup
    )

    consensus = (
        consensus
        .sort_values(
            [
                "consensus_score",
                "structure_coverage_fraction",
                "structure_rank_stability",
                "mean_structure_top_5_fraction",
            ],
            ascending=[
                False,
                False,
                False,
                False,
            ],
        )
        .reset_index(
            drop=True
        )
    )

    consensus.insert(
        0,
        "consensus_rank",
        np.arange(
            1,
            len(
                consensus
            )
            + 1,
        ),
    )

    return (
        consensus,
        structure_scores,
        working,
    )


def select_consensus_visualization_poses(
    dock_scores: pd.DataFrame,
    consensus: pd.DataFrame,
    top_n_compounds: int,
    minimum_dock_percentile: float,
) -> pd.DataFrame:
    """
    Select one representative pose for each top consensus compound.

    Prefer the best GlideScore among docks where the compound
    achieved the configured within-dock percentile. If no such
    dock exists, use the highest-percentile result.
    """
    selected_rows = []

    top_compounds = (
        consensus
        .head(
            top_n_compounds
        )[
            [
                "compound_id",
                "consensus_rank",
                "consensus_score",
                "structure_coverage_fraction",
                "structure_rank_stability",
            ]
        ]
    )

    for consensus_row in (
        top_compounds.itertuples(
            index=False
        )
    ):
        compound_rows = dock_scores[
            dock_scores[
                "compound_id"
            ].astype(
                str
            )
            == str(
                consensus_row.compound_id
            )
        ].copy()

        preferred = compound_rows[
            compound_rows[
                "dock_percentile_score"
            ]
            >= minimum_dock_percentile
        ]

        if preferred.empty:
            selected = (
                compound_rows
                .sort_values(
                    [
                        "dock_percentile_score",
                        "glide_score",
                    ],
                    ascending=[
                        False,
                        True,
                    ],
                )
                .iloc[
                    0
                ]
            )

        else:
            selected = (
                preferred
                .sort_values(
                    [
                        "glide_score",
                        "dock_percentile_score",
                    ],
                    ascending=[
                        True,
                        False,
                    ],
                )
                .iloc[
                    0
                ]
            )

        selected = selected.to_dict()

        selected.update(
            {
                "consensus_rank": int(
                    consensus_row.consensus_rank
                ),
                "consensus_score": float(
                    consensus_row.consensus_score
                ),
                "structure_coverage_fraction": float(
                    consensus_row.structure_coverage_fraction
                ),
                "structure_rank_stability": float(
                    consensus_row.structure_rank_stability
                ),
            }
        )

        selected_rows.append(
            selected
        )

    return pd.DataFrame(
        selected_rows
    )


def save_consensus_score_coverage_plot(
    consensus: pd.DataFrame,
    output_file: Path,
    top_n_labels: int,
):
    figure, axis = plt.subplots(
        figsize=(
            9.0,
            7.0,
        )
    )

    scatter = axis.scatter(
        consensus[
            "consensus_score"
        ],
        consensus[
            "structure_coverage_fraction"
        ],
        c=consensus[
            "structure_rank_stability"
        ],
        s=65,
        cmap="viridis",
        vmin=0.0,
        vmax=1.0,
        edgecolor="black",
        linewidth=0.4,
        alpha=0.85,
    )

    for row in (
        consensus
        .head(
            top_n_labels
        )
        .itertuples(
            index=False
        )
    ):
        axis.annotate(
            str(
                row.compound_id
            ),
            (
                row.consensus_score,
                row.structure_coverage_fraction,
            ),
            xytext=(
                4,
                4,
            ),
            textcoords="offset points",
            fontsize=8,
        )

    axis.set_xlabel(
        "Structure-balanced MCSS consensus score"
    )

    axis.set_ylabel(
        "Structure coverage fraction"
    )

    axis.set_xlim(
        0.0,
        1.02,
    )

    axis.set_ylim(
        0.0,
        1.05,
    )

    axis.set_title(
        "MCSS consensus rank strength and coverage"
    )

    colour_bar = figure.colorbar(
        scatter,
        ax=axis,
    )

    colour_bar.set_label(
        "Rank stability across structures"
    )

    figure.tight_layout()

    figure.savefig(
        output_file,
        dpi=250,
        bbox_inches="tight",
    )

    plt.close(
        figure
    )


def save_top5_recurrence_plot(
    consensus: pd.DataFrame,
    output_file: Path,
    top_n: int,
):
    display = (
        consensus
        .head(
            top_n
        )
        .sort_values(
            "mean_structure_top_5_fraction",
            ascending=True,
        )
    )

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
            10.0,
            figure_height,
        )
    )

    axis.barh(
        np.arange(
            len(
                display
            )
        ),
        display[
            "mean_structure_top_5_fraction"
        ],
        color="#0077B6",
        alpha=0.85,
    )

    axis.set_yticks(
        np.arange(
            len(
                display
            )
        )
    )

    axis.set_yticklabels(
        display[
            "compound_id"
        ]
    )

    axis.set_xlim(
        0.0,
        1.0,
    )

    axis.set_xlabel(
        "Mean within-structure fraction of docks ranked top 5"
    )

    axis.set_ylabel(
        "Compound"
    )

    axis.set_title(
        "Repeated top-five MCSS performance"
    )

    figure.tight_layout()

    figure.savefig(
        output_file,
        dpi=250,
        bbox_inches="tight",
    )

    plt.close(
        figure
    )


def save_compound_structure_heatmap(
    structure_scores: pd.DataFrame,
    consensus: pd.DataFrame,
    output_file: Path,
    top_n: int,
):
    compound_order = (
        consensus
        .head(
            top_n
        )[
            "compound_id"
        ]
        .tolist()
    )

    matrix = (
        structure_scores
        .pivot(
            index="compound_id",
            columns="source_name",
            values="structure_median_percentile",
        )
        .reindex(
            compound_order
        )
    )

    figure_width = max(
        8.0,
        0.8
        * len(
            matrix.columns
        )
        + 3.0,
    )

    figure_height = max(
        6.0,
        0.35
        * len(
            matrix.index
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
        matrix.to_numpy(
            dtype=float
        ),
        cmap="viridis",
        vmin=0.0,
        vmax=1.0,
        aspect="auto",
    )

    axis.set_xticks(
        np.arange(
            len(
                matrix.columns
            )
        )
    )

    axis.set_xticklabels(
        matrix.columns,
        rotation=90,
    )

    axis.set_yticks(
        np.arange(
            len(
                matrix.index
            )
        )
    )

    axis.set_yticklabels(
        matrix.index
    )

    axis.set_xlabel(
        "Cryo-EM source structure"
    )

    axis.set_ylabel(
        "Compound"
    )

    axis.set_title(
        "Median within-structure MCSS percentile"
    )

    colour_bar = figure.colorbar(
        image,
        ax=axis,
    )

    colour_bar.set_label(
        "Median percentile: 1 = best"
    )

    figure.tight_layout()

    figure.savefig(
        output_file,
        dpi=250,
        bbox_inches="tight",
    )

    plt.close(
        figure
    )


def analyse_mcss_consensus(
    complex_jobs,
    poseviewers,
    output_dir: Path,
    config: dict,
):
    """
    Complete activity-independent MCSS consensus analysis.
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

    top_ks = tuple(
        int(
            value
        )
        for value in config.get(
            "top_ks",
            [
                1,
                5,
                10,
            ],
        )
    )

    top_n_plot = int(
        config.get(
            "top_n_plot",
            30,
        )
    )

    top_n_labels = int(
        config.get(
            "top_n_labels",
            20,
        )
    )

    top_n_pymol = int(
        config.get(
            "top_n_pymol",
            20,
        )
    )

    minimum_dock_percentile = float(
        config.get(
            "visualization_minimum_dock_percentile",
            0.90,
        )
    )

    binding_site_cutoff = float(
        config.get(
            "binding_site_cutoff",
            5.0,
        )
    )

    dock_input = collect_mcss_docking_scores(
        complex_jobs=complex_jobs,
        poseviewers=poseviewers,
    )

    (
        consensus,
        structure_scores,
        dock_scores,
    ) = calculate_mcss_consensus_ranking(
        dock_scores=dock_input,
        top_ks=top_ks,
    )

    dock_scores.to_csv(
        output_dir
        / "mcss_compound_dock_ranks.csv",
        index=False,
    )

    structure_scores.to_csv(
        output_dir
        / "mcss_compound_structure_ranks.csv",
        index=False,
    )

    consensus.to_csv(
        output_dir
        / "mcss_compound_consensus_ranking.csv",
        index=False,
    )

    save_consensus_score_coverage_plot(
        consensus=consensus,
        output_file=(
            output_dir
            / "mcss_consensus_score_vs_coverage.png"
        ),
        top_n_labels=top_n_labels,
    )

    save_top5_recurrence_plot(
        consensus=consensus,
        output_file=(
            output_dir
            / "mcss_consensus_top5_recurrence.png"
        ),
        top_n=top_n_plot,
    )

    save_compound_structure_heatmap(
        structure_scores=structure_scores,
        consensus=consensus,
        output_file=(
            output_dir
            / "mcss_consensus_structure_heatmap.png"
        ),
        top_n=top_n_plot,
    )

    selected_poses = (
        select_consensus_visualization_poses(
            dock_scores=dock_scores,
            consensus=consensus,
            top_n_compounds=top_n_pymol,
            minimum_dock_percentile=(
                minimum_dock_percentile
            ),
        )
    )

    selected_poses.to_csv(
        output_dir
        / "mcss_consensus_selected_poses.csv",
        index=False,
    )

    pymol_root = (
        output_dir
        / "consensus_pymol"
    )

    pymol_root.mkdir(
        parents=True,
        exist_ok=True,
    )

    pymol_scripts = []

    for row in selected_poses.itertuples(
        index=False
    ):
        try:
            script = export_consensus_compound_pose(
                poseviewer=Path(
                    row.poseviewer_path
                ),
                compound_id_value=str(
                    row.compound_id
                ),
                compound_pose_index=int(
                    row.compound_pose_index
                ),
                reference_ligand=Path(
                    row.reference_ligand_path
                ),
                output_dir=(
                    pymol_root
                    / (
                        f"rank{int(row.consensus_rank):03d}_"
                        f"{str(row.compound_id)}"
                    )
                ),

                binding_site_cutoff=(
                    binding_site_cutoff
                ),

                metadata={
                    "consensus_rank": int(
                        row.consensus_rank
                    ),
                    "consensus_score": float(
                        row.consensus_score
                    ),
                    "structure_coverage_fraction": float(
                        row.structure_coverage_fraction
                    ),
                    "structure_rank_stability": float(
                        row.structure_rank_stability
                    ),
                    "source_name": str(
                        row.source_name
                    ),
                    "dock_name": str(
                        row.dock_name
                    ),
                    "dock_rank": int(
                        row.dock_rank
                    ),
                    "dock_percentile_score": float(
                        row.dock_percentile_score
                    ),
                    "glide_score": float(
                        row.glide_score
                    ),
                },
            )

            pymol_scripts.append(
                script
            )

        except Exception as exc:
            logger.warning(
                "[%s] Consensus PyMOL export failed: %s",
                row.compound_id,
                exc,
                exc_info=True,
            )

    # --------------------------------------------------
    # Combined PyMOL bundle
    #
    # One selected consensus ligand + its associated receptor
    # per ranked compound. No reference ligand is included.
    # --------------------------------------------------

    combined_pymol_dir = (
        output_dir
        / "consensus_pymol_combined"
    )

    combined_pymol_script = None

    try:
        combined_pymol_script = (
            export_combined_consensus_poses(
                selected_poses=selected_poses,
                output_dir=(
                    combined_pymol_dir
                ),
                binding_site_cutoff=(
                    binding_site_cutoff
                ),
            )
        )

    except Exception as exc:
        logger.warning(
            "Combined consensus PyMOL export failed: %s",
            exc,
            exc_info=True,
        )

    else:
        logger.info(
            "Combined consensus PyMOL bundle generated: %s",
            combined_pymol_script,
        )

    logger.info(
        "Generated MCSS consensus ranking for %d compounds "
        "across %d structures and %d MCSS docks",
        len(
            consensus
        ),
        structure_scores[
            "source_name"
        ].nunique(),
        dock_scores[
            "dock_name"
        ].nunique(),
    )

    return {
        "consensus_ranking": (
            output_dir
            / "mcss_compound_consensus_ranking.csv"
        ),
        "structure_ranking": (
            output_dir
            / "mcss_compound_structure_ranks.csv"
        ),
        "dock_ranking": (
            output_dir
            / "mcss_compound_dock_ranks.csv"
        ),
        "selected_poses": (
            output_dir
            / "mcss_consensus_selected_poses.csv"
        ),
        "pymol_scripts": pymol_scripts,

        "combined_pymol_script": (
            combined_pymol_script
        ),
    }