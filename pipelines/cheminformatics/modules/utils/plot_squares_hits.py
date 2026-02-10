import pandas as pd
import numpy as np
import seaborn as sns

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
# from matplotlib.gridspec import GridSpec

from pathlib import Path
# from upsetplot import from_indicators, UpSet
# from scipy.stats import gaussian_kde
from itertools import combinations

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def plot_squares_hits(df, top_n=20, enrichment_col="enrichment_mean",
                          enrichment_threshold=1.5, thresholds=None, output_dir=None):
    """
    Comprehensive DEL squares plotting including barplots, heatmaps, scatter, and threshold sweep.
    
    Parameters
    ----------
    df : pd.DataFrame
        Squares hits with columns ['NsynthonID', 'condition', enrichment_col, 'squares_score'].
    top_n : int
        Number of top synthons to show in barplots.
    enrichment_col : str
        Column to use for enrichment.
    enrichment_threshold : float
        Threshold to define “active” for squares.
    thresholds : list[float]
        Thresholds to sweep for visualization of effect on hits/squares.
    output_dir : str | Path
        Directory to save plots.
    """


    output_dir = Path(output_dir) if output_dir is not None else Path(".")
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_files = {}

    # ----------------------------
    # Top synthons by squares_score
    # ----------------------------
    top_synthons = df.groupby("NsynthonID", observed=True)["squares_score"].max().nlargest(top_n).index
    df_top = df[df["NsynthonID"].isin(top_synthons)].copy()
    syn_order = df_top.groupby("NsynthonID", observed=True)["squares_score"].max().sort_values(ascending=False).index
    df_top["NsynthonID"] = pd.Categorical(df_top["NsynthonID"], categories=syn_order, ordered=True)
    mean_enrichment = df_top.groupby("NsynthonID", observed=True)[enrichment_col].mean()

    def colored_barplot(values, categories, mean_enr, cmap_name, ylabel, title, filename):
        norm = mcolors.Normalize(vmin=mean_enr.min(), vmax=mean_enr.max())
        cmap = cm.get_cmap(cmap_name)
        colors = cmap(norm(mean_enr.loc[categories]))
        fig, ax = plt.subplots(figsize=(12,6))
        bars = ax.bar(categories, values, color=colors)
        ax.set_xticklabels(categories, rotation=45, ha="right")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, label="Mean Enrichment")
        plt.tight_layout()
        out_path = output_dir / filename
        fig.savefig(out_path, dpi=300)
        plt.close()
        return out_path

    # Total squares
    total_squares = df_top.groupby("NsynthonID", observed=True)["squares_score"].max()
    plot_files["total_squares"] = colored_barplot(
        total_squares.values, total_squares.index, mean_enrichment,
        cmap_name="viridis",
        ylabel="Total Squares Score",
        title="Top Synthons by Total Squares Score (colored by mean enrichment)",
        filename="squares_total_score_barplot.png"
    )

    # Max enrichment
    max_enr = df_top.groupby("NsynthonID", observed=True)[enrichment_col].max()
    plot_files["max_enrichment"] = colored_barplot(
        max_enr.values, max_enr.index, mean_enrichment,
        cmap_name="magma",
        ylabel="Max Enrichment",
        title="Max Condition Enrichment per Synthon (colored by mean enrichment)",
        filename="squares_max_enrichment_barplot.png"
    )

    # Active conditions
    active_conditions = df_top.groupby("NsynthonID", observed=True).apply(
        lambda g: (g[enrichment_col] >= enrichment_threshold).sum()
    )
    plot_files["active_conditions"] = colored_barplot(
        active_conditions.values, active_conditions.index, mean_enrichment,
        cmap_name="cividis",
        ylabel=f"Number of Active Conditions (≥ {enrichment_threshold})",
        title="Breadth of Activity per Synthon (colored by mean enrichment)",
        filename="squares_active_conditions_barplot.png"
    )

    # Mean std enrichment
    mean_std_enr = df_top.groupby("NsynthonID", observed=True)[enrichment_col].agg(["mean","std"])
    fig, ax = plt.subplots(figsize=(12,6))
    ax.bar(mean_std_enr.index, mean_std_enr["mean"].values, yerr=mean_std_enr["std"].values, color="skyblue", capsize=4)
    ax.set_xticklabels(mean_std_enr.index, rotation=45, ha="right")
    ax.set_ylabel("Mean Enrichment ± Std")
    ax.set_title("Mean and Variability of Enrichment per Synthon")
    plt.tight_layout()
    plot_files["mean_std_enrichment"] = output_dir / "squares_mean_std_enrichment_barplot.png"
    fig.savefig(plot_files["mean_std_enrichment"], dpi=300)
    plt.close()

    # Squares contribution heatmap
    conditions = df_top["condition"].unique()
    contribution_matrix = pd.DataFrame(0, index=df_top["NsynthonID"].cat.categories, columns=conditions)
    for syn in df_top["NsynthonID"].cat.categories:
        syn_df = df_top[df_top["NsynthonID"] == syn]
        active_conditions_syn = syn_df[syn_df[enrichment_col] >= enrichment_threshold]["condition"].tolist()
        for c1, c2 in combinations(active_conditions_syn, 2):
            contribution_matrix.loc[syn, c1] += 1
            contribution_matrix.loc[syn, c2] += 1
    fig, ax = plt.subplots(figsize=(12, max(6, len(contribution_matrix)*0.25)))
    sns.heatmap(contribution_matrix, cmap="viridis", linewidths=0.5, linecolor="gray", ax=ax)
    ax.set_ylabel("Synthon")
    ax.set_xlabel("Condition")
    ax.set_title("Per-Condition Contribution to Squares Score")
    plt.tight_layout()
    plot_files["squares_contribution_heatmap"] = output_dir / "squares_contribution_heatmap.png"
    fig.savefig(plot_files["squares_contribution_heatmap"], dpi=300)
    plt.close()

    # Scatter squares_score vs mean enrichment
    fig, ax = plt.subplots(figsize=(8,6))
    sns.scatterplot(x=total_squares.values, y=mean_enrichment.values, ax=ax)
    for i, syn in enumerate(total_squares.index):
        ax.text(total_squares.values[i], mean_enrichment.values[i], syn, fontsize=8)
    ax.set_xlabel("Squares Score")
    ax.set_ylabel("Mean Enrichment")
    ax.set_title("Squares Score vs Mean Enrichment")
    plt.tight_layout()
    plot_files["scatter_score_vs_mean"] = output_dir / "squares_score_vs_mean_enrichment.png"
    fig.savefig(plot_files["scatter_score_vs_mean"], dpi=300)
    plt.close()

    # Binary heatmap: active/inactive
    binary_df = df_top.copy()
    binary_df["active"] = binary_df[enrichment_col] >= enrichment_threshold
    binary_matrix = binary_df.pivot(index="NsynthonID", columns="condition", values="active").fillna(False)
    fig, ax = plt.subplots(figsize=(12, max(6, len(binary_matrix)*0.25)))
    sns.heatmap(binary_matrix.astype(int), cmap="Greys", linewidths=0.5, linecolor="gray", cbar=False, ax=ax)
    ax.set_ylabel("Synthon")
    ax.set_xlabel("Condition")
    ax.set_title(f"Binary Active/Inactive per Condition (Threshold ≥ {enrichment_threshold})")
    plt.tight_layout()
    plot_files["binary_heatmap"] = output_dir / "squares_binary_heatmap.png"
    fig.savefig(plot_files["binary_heatmap"], dpi=300)
    plt.close()

    # ----------------------------
    # Threshold sweep plots
    # ----------------------------
    if thresholds is None:
        thresholds = np.linspace(1.0, 3.0, 11)

    active_counts = []
    total_squares_list = []
    for th in thresholds:
        df["active"] = df[enrichment_col] >= th
        active_counts.append(df.groupby("NsynthonID", observed=True)["active"].any().sum())
        squares_list = []
        conds = df["condition"].unique()
        for c1, c2 in combinations(conds, 2):
            df_c1 = df[df["condition"] == c1][["NsynthonID","active"]].set_index("NsynthonID")
            df_c2 = df[df["condition"] == c2][["NsynthonID","active"]].set_index("NsynthonID")
            squares_list.extend(set(df_c1[df_c1["active"]].index) & set(df_c2[df_c2["active"]].index))
        total_squares_list.append(len(squares_list))

    # Active synthons vs threshold
    fig, ax = plt.subplots(figsize=(8,5))
    ax.plot(thresholds, active_counts, marker='o', label="Active Synthons")
    ax.set_xlabel("Enrichment Threshold")
    ax.set_ylabel("Number of Active Synthons")
    ax.set_title("Active Synthons vs Enrichment Threshold")
    ax.grid(True)
    ax.legend()
    plt.tight_layout()
    plot_files["threshold_sweep_active_synthons"] = output_dir / "threshold_sweep_active_synthons.png"
    fig.savefig(plot_files["threshold_sweep_active_synthons"], dpi=300)
    plt.close()

    # Total squares vs threshold
    fig, ax = plt.subplots(figsize=(8,5))
    ax.plot(thresholds, total_squares_list, marker='o', color='orange', label="Total Squares")
    ax.set_xlabel("Enrichment Threshold")
    ax.set_ylabel("Total Squares")
    ax.set_title("Total Squares vs Enrichment Threshold")
    ax.grid(True)
    ax.legend()
    plt.tight_layout()
    plot_files["threshold_sweep_total_squares"] = output_dir / "threshold_sweep_total_squares.png"
    fig.savefig(plot_files["threshold_sweep_total_squares"], dpi=300)
    plt.close()

    return plot_files




