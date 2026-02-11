import pandas as pd
import numpy as np
import seaborn as sns

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec

from pathlib import Path
# from upsetplot import from_indicators, UpSet
# from scipy.stats import gaussian_kde
from itertools import combinations

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def plot_squares_hits(df, full_df=None, top_n=20, enrichment_col="enrichment_mean",
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
        ax.set_xticks(range(len(categories)))
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

    # ----------------------------
    # Normalised enrichment heatmap (Z-scored per condition)
    # ----------------------------

    # Decide which full dataframe to use for enrichment values
    enr_df = full_df if full_df is not None else df

    # Keep only synthons in df_top
    plot_syn_ids = df_top["NsynthonID"].unique()
    heatmap_df = enr_df[enr_df["NsynthonID"].isin(plot_syn_ids)].copy()

    # Pivot to matrix: rows = synthon, cols = condition
    pivot_df = heatmap_df.pivot(index="NsynthonID", columns="condition", values=enrichment_col)

    # Compute per-condition z-score using all synthons
    pivot_z = pivot_df.copy()
    for col in pivot_z.columns:
        col_values = enr_df[enr_df["condition"] == col][enrichment_col]
        col_mean = col_values.mean()
        col_std = col_values.std()
        if col_std > 0:
            pivot_z[col] = (pivot_z[col] - col_mean) / col_std
        else:
            pivot_z[col] = 0.0

    # Order rows by squares_score from df_top
    row_order = (
        df_top.groupby("NsynthonID", observed=True)["squares_score"]
        .max()
        .sort_values(ascending=False)
        .index
    )
    pivot_z = pivot_z.loc[row_order]

    fig, ax = plt.subplots(figsize=(12, max(6, len(pivot_z)*0.3)))
    sns.heatmap(
        pivot_z,
        cmap="coolwarm",
        center=0,
        linewidths=0.5,
        linecolor="gray",
        ax=ax
    )
    ax.set_ylabel("Synthon")
    ax.set_xlabel("Condition")
    ax.set_title("Top Synthons — Condition-Normalized Enrichment (Z-score)")
    plt.tight_layout()
    plot_files["normalized_enrichment_heatmap"] = output_dir / "squares_normalized_enrichment_heatmap.png"
    fig.savefig(plot_files["normalized_enrichment_heatmap"], dpi=300)
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

    # ----------------------------
    # Binary heatmap: active/inactive (full conditions, top synthons only)
    # ----------------------------
    # Use full_df if available
    bin_df = full_df if full_df is not None else df

    # Keep only synthons in df_top
    bin_df = bin_df[bin_df["NsynthonID"].isin(df_top["NsynthonID"].unique())].copy()

    # Compute active/inactive based on enrichment threshold
    bin_df["active"] = bin_df[enrichment_col] >= enrichment_threshold

    # Pivot to matrix: rows = synthon, cols = condition
    binary_matrix = bin_df.pivot(index="NsynthonID", columns="condition", values="active").fillna(False)

    # Order rows to match normalized heatmap
    binary_matrix = binary_matrix.loc[row_order]

    inactive_matrix = (~binary_matrix).astype(int)

    fig, ax = plt.subplots(figsize=(12, max(6, len(binary_matrix)*0.25)))
    sns.heatmap(
        inactive_matrix.astype(int),
        cmap="Blues",        # diverging or sequential colormap
        linewidths=0.5,
        linecolor="gray",
        cbar=False,
        ax=ax,
        alpha=0.7
    )
    ax.set_ylabel("Synthon")
    ax.set_xlabel("Condition")
    ax.set_title(f"Binary Active/Inactive per Condition (Threshold ≥ {enrichment_threshold})")
    plt.tight_layout()
    plot_files["binary_heatmap"] = output_dir / "squares_binary_heatmap.png"
    fig.savefig(plot_files["binary_heatmap"], dpi=300)
    plt.close()

    # ----------------------------
    # Threshold sweep plots: How many hits do we lose with different activity thresholds?
    # ----------------------------
    if thresholds is None:
        thresholds = np.linspace(1.0, 3.0, 11)

    active_counts = []
    total_squares_list = []

    sweep_df = full_df if full_df is not None else df

    for th in thresholds:
        tmp = sweep_df.copy()
        tmp["active"] = tmp[enrichment_col] >= th

        active_counts.append(
            tmp.groupby("NsynthonID", observed=True)["active"].any().sum()
        )

        squares_list = []
        conds = tmp["condition"].unique()

        for c1, c2 in combinations(conds, 2):
            df_c1 = tmp[tmp["condition"] == c1][["NsynthonID", "active"]].set_index("NsynthonID")
            df_c2 = tmp[tmp["condition"] == c2][["NsynthonID", "active"]].set_index("NsynthonID")
            squares_list.extend(
                set(df_c1[df_c1["active"]].index)
                & set(df_c2[df_c2["active"]].index)
            )

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

    # ----------------------------
    # Faceted UMAP by condition (highlight top squares hits)
    # ----------------------------
    if {"umap_1", "umap_2"}.issubset(df.columns):
        # Optional: library-weighted subsampling for plotting
        df_for_umap = df.copy()
        if "library" in df_for_umap.columns:
            lib_counts = df_for_umap["library"].value_counts()
            df_for_umap["umap_weight"] = 1.0 / lib_counts[df_for_umap["library"]].values
            # Sample all rows proportionally to library weight
            df_for_umap = df_for_umap.sample(n=len(df_for_umap), weights="umap_weight", random_state=42).reset_index(drop=True)

        plot_col = "squares_score"
        conditions = sorted(df_for_umap["condition"].unique())
        n_cols = 3
        n_rows = int(np.ceil(len(conditions) / n_cols))

        fig = plt.figure(figsize=(4 * n_cols, 4 * n_rows + 1))
        gs = GridSpec(n_rows + 1, n_cols, height_ratios=[4]*n_rows + [0.2], hspace=0.3)

        vmin = df_for_umap[plot_col].min()
        vmax = df_for_umap[plot_col].max()
        for i, cond in enumerate(conditions):
            row = i // n_cols
            col = i % n_cols
            ax = fig.add_subplot(gs[row, col])

            # background points
            ax.scatter(df_for_umap["umap_1"], df_for_umap["umap_2"], s=10, c="lightgray", alpha=0.2, linewidth=0)

            # highlight points for this condition
            sub = df_for_umap[df_for_umap["condition"] == cond]
            sc = ax.scatter(sub["umap_1"], sub["umap_2"], c=sub[plot_col], cmap="plasma_r",
                            vmin=vmin, vmax=vmax, s=20, alpha=0.9, linewidth=0)

            ax.set_title(cond)
            ax.set_xticks([])
            ax.set_yticks([])

        for i in range(len(conditions), n_rows*n_cols):
            ax = fig.add_subplot(gs[i // n_cols, i % n_cols])
            ax.axis("off")

        # horizontal colorbar
        cbar_ax = fig.add_subplot(gs[-1, :])
        plt.colorbar(sc, cax=cbar_ax, orientation="horizontal", label=plot_col)
        fig.suptitle(f"UMAP of DEL hits (highlighted by {plot_col})", y=1.02)
        umap_file = Path(output_dir) / "umap_by_condition.png"
        plt.savefig(umap_file, dpi=300, bbox_inches="tight")
        plt.close()
        plot_files["umap_by_condition"] = umap_file
    else:
        logger.warning("[plot_squares_hits] UMAP scatter skipped: umap_1 / umap_2 not found in dataframe")


    return plot_files




