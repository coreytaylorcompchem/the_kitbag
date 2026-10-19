import math

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.colors import (
    BoundaryNorm,
    ListedColormap,
)
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import dendrogram

from scipy.stats import kendalltau, pearsonr, spearmanr

from schrodinger import structure

from modules.utils.glideem_helpers import (
    glide_smiles_by_lignum,
    glide_score,
    glide_score_components,
    property_value,
    compound_id,
    smiles_value,
)

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def analyse_correlations(
    complex_jobs,
    poseviewers,
    ligand_file,
    source_ligand_file,
    activity_column,
    output_dir,
    min_pairs=5,
):
    activities = {}
    ligand_metadata = {}

    for st in structure.StructureReader(
        str(ligand_file)
    ):
        cid = compound_id(st)

        if not cid:
            continue

        try:
            value = float(
                property_value(
                    st,
                    activity_column,
                )
            )
        except (TypeError, ValueError):
            continue

        if not math.isfinite(value):
            continue

        activities.setdefault(
            cid,
            value,
        )

        ligand_metadata.setdefault(
            cid,
            {
                "smiles": smiles_value(st),
            },
        )
    if not activities:
        raise RuntimeError(f"No numeric values found for activity property '{activity_column}'")

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

        smiles_by_lignum = glide_smiles_by_lignum(
            pv
        )

        for st in structure.StructureReader(
            str(pv)
        ):
            cid = compound_id(st)
            score = glide_score(st)

            lignum = st.property.get(
                "i_i_glide_lignum"
            )

            try:
                lignum = int(lignum)
            except (TypeError, ValueError):
                lignum = None

            if cid:
                docked_ids.add(cid)

            if cid and score is not None:
                scored_ids.add(cid)

            if (
                cid in activities
                and score is not None
                and (
                    cid not in best
                    or score < best[cid]["score"]
                )
            ):
                best[cid] = {
                    "score": score,
                    "structure": st,
                    "lignum": lignum,
                }

                matched_ids.add(cid)

        logger.debug(
            "[%s] Docked IDs with scores: %d; "
            "matched to activity: %d",
            name,
            len(best),
            sum(
                cid in activities
                for cid in best
            ),
        )

        for cid, result in best.items():

            score = result["score"]
            docked_st = result["structure"]

            lignum = result["lignum"]

            smiles = (
                smiles_by_lignum.get(lignum)
                if lignum is not None
                else None
            )

            components = glide_score_components(
                docked_st
            )

            metadata = ligand_metadata.get(
                cid,
                {}
            )

            details.append(
                {
                    "complex_pose": name,
                    "source_name": by_name[name]["source_name"],
                    "pose_number": by_name[name]["pose_number"],

                    "compound_id": cid,
                    "smiles": smiles,

                    "experimental_potency": activities[cid],
                    "activity_column": activity_column,

                    "glide_score": score,

                    **components,

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
    
    # if not details:
    #     raise RuntimeError("No matched docking score/activity pairs")

    output_dir.mkdir(parents=True, exist_ok=True)

    detail_df = pd.DataFrame(details)
    detail_df.to_csv(output_dir / "matched_scores_and_activities.csv", index=False)
    summaries = []
    for name, group in detail_df.groupby("complex_pose"):
        if len(group) < min_pairs:
            logger.warning("%s has only %d matched pairs; excluded", name, len(group))
            continue
        x = group["glide_score"].to_numpy(float)
        y = group["experimental_potency"].to_numpy(float)
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

        source_name = group[
                "source_name"
            ].iloc[0]

        pose_number = int(
            group[
                "pose_number"
            ].iloc[0]
        )

        plot_dir = (
            output_dir /
            source_name
        )

        plot_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        pose_table = group[
            [
                "compound_id",
                "smiles",
                "experimental_potency",
                "glide_score",
                "docking_score",
                "glide_emodel",
                "glide_energy",
                "glide_evdw",
                "glide_ecoul",
                "glide_lipo",
                "glide_hbond",
                "glide_erotb",
                "glide_esite",
                "glide_einternal",
            ]
        ].copy()

        pose_table = pose_table.sort_values(
            "experimental_potency",
            ascending=False,
        )

        pose_table.to_csv(
            plot_dir /
            (
                f"pose{pose_number:02d}"
                "_docking_vs_affinity.csv"
            ),
            index=False,
        )

        # Conventional experimental vs predicted orientation:
        # experimental values on X, predicted values on Y.
        x = group[
            "experimental_potency"
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
            plot_dir /
            (
                f"pose{pose_number:02d}"
                "_experimental_vs_predicted.png"
            ),
            dpi=200,
        )

        plt.close(fig)
    best = summary.iloc[0]
    (output_dir / "best_binding_pose.txt").write_text(
        f"Best complex pose: {best.complex_pose}\nN: {best.n}\nPearson r: {best.pearson_r:.6f}\n"
        f"Spearman rho: {best.spearman_rho:.6f}\nKendall tau: {best.kendall_tau:.6f}\n"
    )
    return summary

def _ensure_parent(
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

    return output_file


def save_clustered_score_vs_rmsd(
    metrics: pd.DataFrame,
    score_column: str,
    output_file: Path,
    title: str,
    y_label: str,
    rmsd_column: str = "Ligand RMSD",
    cluster_column: str = "cluster_id",
    dpi: int = 200,
) -> Path:
    """
    Plot one GlideEM score against RMSD, coloured by cluster.
    """

    output_file = _ensure_parent(
        output_file
    )

    required = {
        rmsd_column,
        score_column,
        cluster_column,
    }

    missing = required - set(
        metrics.columns
    )

    if missing:
        raise ValueError(
            "Missing plotting columns: "
            f"{sorted(missing)}"
        )

    plot_data = (
        metrics[
            [
                rmsd_column,
                score_column,
                cluster_column,
            ]
        ]
        .dropna()
        .copy()
    )

    if plot_data.empty:
        raise ValueError(
            f"No finite data available for {title}"
        )

    cluster_ids = sorted(
        plot_data[
            cluster_column
        ]
        .astype(int)
        .unique()
    )

    cluster_to_index = {
        cluster_id: index
        for index, cluster_id in enumerate(
            cluster_ids
        )
    }

    colour_indices = (
        plot_data[
            cluster_column
        ]
        .astype(int)
        .map(cluster_to_index)
        .to_numpy(int)
    )

    n_clusters = len(
        cluster_ids
    )

    cmap = plt.get_cmap(
        "tab20",
        max(n_clusters, 1),
    )

    norm = BoundaryNorm(
        np.arange(
            -0.5,
            n_clusters + 0.5,
            1.0,
        ),
        cmap.N,
    )

    fig, ax = plt.subplots(
        figsize=(7, 6)
    )

    scatter = ax.scatter(
        plot_data[
            rmsd_column
        ],
        plot_data[
            score_column
        ],
        c=colour_indices,
        cmap=cmap,
        norm=norm,
        alpha=0.85,
    )

    ax.set_xlabel(
        "RMSD from input cryo-EM ligand pose (Å)"
    )

    ax.set_ylabel(
        y_label
    )

    ax.set_title(
        title
    )

    colour_bar = fig.colorbar(
        scatter,
        ax=ax,
        ticks=np.arange(
            n_clusters
        ),
    )

    colour_bar.set_label(
        "Pose cluster"
    )

    colour_bar.ax.set_yticklabels(
        [
            f"C{cluster_id}"
            for cluster_id in cluster_ids
        ]
    )

    fig.tight_layout()

    fig.savefig(
        output_file,
        dpi=dpi,
    )

    plt.close(fig)

    return output_file


def save_rmsd_heatmap(
    rmsd_matrix: np.ndarray,
    pose_labels: list[str],
    output_file: Path,
    title: str,
    dpi: int = 200,
) -> Path:
    output_file = _ensure_parent(
        output_file
    )

    rmsd_matrix = np.asarray(
        rmsd_matrix,
        dtype=float,
    )

    expected_shape = (
        len(pose_labels),
        len(pose_labels),
    )

    if rmsd_matrix.shape != expected_shape:
        raise ValueError(
            "RMSD matrix shape does not match labels: "
            f"{rmsd_matrix.shape} versus "
            f"{expected_shape}"
        )

    fig, ax = plt.subplots(
        figsize=(8, 7)
    )

    image = ax.imshow(
        rmsd_matrix,
        interpolation="nearest",
        aspect="auto",
    )

    positions = np.arange(
        len(pose_labels)
    )

    ax.set_xticks(
        positions
    )

    ax.set_yticks(
        positions
    )

    ax.set_xticklabels(
        pose_labels,
        rotation=90,
    )

    ax.set_yticklabels(
        pose_labels,
    )

    ax.set_title(
        title
    )

    fig.colorbar(
        image,
        ax=ax,
        label="Heavy-atom RMSD (Å)",
    )

    fig.tight_layout()

    fig.savefig(
        output_file,
        dpi=dpi,
    )

    plt.close(fig)

    return output_file


def save_cluster_dendrogram(
    linkage_matrix: np.ndarray,
    pose_labels: list[str],
    cluster_cutoff: float,
    output_file: Path,
    title: str,
    dpi: int = 200,
) -> Path:
    output_file = _ensure_parent(
        output_file
    )

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
        title
    )

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        output_file,
        dpi=dpi,
    )

    plt.close(fig)

    return output_file


def save_pose_ambiguity_map(
    summary_df: pd.DataFrame,
    output_file: Path,
    dpi: int = 200,
) -> Path:
    output_file = _ensure_parent(
        output_file
    )

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
        output_file,
        dpi=dpi,
    )

    plt.close(fig)

    return output_file


def save_reference_rmsd_distributions(
    all_metrics_df: pd.DataFrame,
    job_order: list[str],
    output_file: Path,
    dpi: int = 200,
) -> Path:
    output_file = _ensure_parent(
        output_file
    )

    distributions = [
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
        distributions,
        positions=positions,
        vert=False,
        showmeans=False,
        showmedians=True,
        showextrema=True,
    )

    for body in violin_parts["bodies"]:
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
        output_file,
        dpi=dpi,
    )

    plt.close(fig)

    return output_file


def save_cluster_composition(
    composition_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    job_order: list[str],
    output_file: Path,
    minimum_label_fraction: float = 0.07,
    dpi: int = 200,
) -> Path:
    output_file = _ensure_parent(
        output_file
    )

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
        row["source_name"]: int(
            row["n_clusters"]
        )
        for _, row in summary_df.iterrows()
    }

    ylabels = [
        (
            f"{job_name}  "
            f"(k={clusters_per_job[job_name]})"
        )
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
        values = (
            composition_pivot[
                cluster_id
            ].to_numpy(float)
        )

        bars = ax.barh(
            composition_pivot.index,
            values,
            left=left,
        )

        for bar, value, start in zip(
            bars,
            values,
            left,
        ):
            if value >= minimum_label_fraction:
                ax.text(
                    start + value / 2,
                    (
                        bar.get_y()
                        + bar.get_height() / 2
                    ),
                    f"C{cluster_id}",
                    ha="center",
                    va="center",
                    fontsize=8,
                )

        left += values

    positions = np.arange(
        len(composition_pivot)
    )

    ax.set_yticks(
        positions
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
        output_file,
        dpi=dpi,
    )

    plt.close(fig)

    return output_file


def save_density_rmsd_small_multiples(
    all_metrics_df: pd.DataFrame,
    job_order: list[str],
    output_file: Path,
    top_n_clusters: int = 10,
    dpi: int = 200,
) -> Path:
    """
    Generate clustered density-score/RMSD plots.

    For each job, the most populous clusters are individually
    coloured. Remaining clusters are grouped as Other.
    """

    output_file = _ensure_parent(
        output_file
    )

    n_jobs = len(
        job_order
    )

    if n_jobs == 0:
        raise ValueError(
            "No GlideEM jobs supplied for plotting."
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
            .dropna(
                subset=[
                    "Ligand RMSD",
                    "Glide Density Score",
                    "cluster_id",
                ]
            )
            .copy()
        )

        if job_metrics.empty:
            axis.set_visible(
                False
            )
            continue

        cluster_counts = (
            job_metrics[
                "cluster_id"
            ]
            .astype(int)
            .value_counts()
        )

        top_cluster_ids = (
            cluster_counts
            .head(top_n_clusters)
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

        if has_other:
            other_index = 0

            cluster_to_colour_index = {
                cluster_id: index
                for index, cluster_id in enumerate(
                    top_cluster_ids,
                    start=1,
                )
            }

        else:
            other_index = None

            cluster_to_colour_index = {
                cluster_id: index
                for index, cluster_id in enumerate(
                    top_cluster_ids
                )
            }

        def colour_index(
            cluster_id,
        ) -> int:
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
            .map(colour_index)
            .to_numpy(int)
        )

        base_cmap = plt.get_cmap(
            "tab10"
        )

        if has_other:
            colours = [
                "#bdbdbd"
            ]

            colours.extend(
                base_cmap(index)
                for index in range(
                    len(top_cluster_ids)
                )
            )

        else:
            colours = [
                base_cmap(index)
                for index in range(
                    len(top_cluster_ids)
                )
            ]

        cmap = ListedColormap(
            colours
        )

        n_categories = len(
            colours
        )

        norm = BoundaryNorm(
            np.arange(
                -0.5,
                n_categories + 0.5,
                1.0,
            ),
            cmap.N,
        )

        axis.scatter(
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

        legend_handles = []

        for cluster_id in top_cluster_ids:
            colour_position = (
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
                        norm(
                            colour_position
                        )
                    ),
                    edgecolor="none",
                    label=(
                        f"C{cluster_id}: "
                        f"{count} poses"
                    ),
                )
            )

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

    for axis in axes.flat[n_jobs:]:
        axis.set_visible(
            False
        )

    fig.suptitle(
        "Clustered GlideEM density-score × RMSD",
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
        output_file,
        dpi=dpi,
        bbox_inches="tight",
    )

    plt.close(fig)

    return output_file