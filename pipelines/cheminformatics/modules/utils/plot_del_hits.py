import pandas as pd
import numpy as np
import seaborn as sns

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

from pathlib import Path
from upsetplot import from_indicators, UpSet
from scipy.stats import gaussian_kde
from itertools import combinations

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def reduce_hits_for_plotting(
    df: pd.DataFrame,
    max_hits: int = 300,
    sort_cols=("p_active", "beta_mean"),
):
    """
    Reduce hits for plotting while preserving cluster diversity.
    """
    if len(df) <= max_hits:
        return df.copy()

    if "cluster_id" not in df.columns:
        return (
            df.sort_values(list(sort_cols), ascending=[False, False])
            .head(max_hits)
            .copy()
        )

    clusters = df["cluster_id"].unique()
    clusters = [c for c in clusters if c != -1]

    if not clusters:
        return (
            df.sort_values(list(sort_cols), ascending=[False, False])
            .head(max_hits)
            .copy()
        )

    n_per = max(1, max_hits // len(clusters))
    pieces = []

    for c in clusters:
        sub = (
            df[df["cluster_id"] == c]
            .sort_values(list(sort_cols), ascending=[False, False])
            .head(n_per)
        )
        pieces.append(sub)

    reduced = (
        pd.concat(pieces)
        .sort_values(list(sort_cols), ascending=[False, False])
        .head(max_hits)
    )

    return reduced.copy()

def plot_del_hits(
    hits_df: pd.DataFrame,
    output_dir: str = "outputs/plots",
    top_n: int = 10,
    physchem_cols=None,
    heatmap_top_n: int = 100,
    full_hits_df: pd.DataFrame | None = None,
):
    """
    Generate DEL hit visualizations:
    1. Heatmap of p_active per synthon x condition (top N synthons by max p_active)
    2. Strip plot of p_active per condition
    3. Bar plot of top hits per condition
    4. UpSet plot of hits overlap across conditions
    5. Scatterplots of p_active vs physchem properties, colored by condition
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_files = {}

    hits_df = hits_df.copy()
    hits_df["condition"] = hits_df["condition"].astype(str)

    if full_hits_df is None:
        full_hits_df = hits_df

    plot_col = "p_active_cond" if "p_active_cond" in hits_df.columns else "p_active"
    hits_df[plot_col] = pd.to_numeric(hits_df[plot_col], errors="coerce")

    if physchem_cols is None:
        physchem_cols = [
            c for c in hits_df.columns
            if c not in {"NsynthonID", "condition", "beta_mean", "beta_hdi_lower", "beta_hdi_upper", "p_active", "p_active_cond"}
            and pd.api.types.is_numeric_dtype(hits_df[c])
        ]
    physchem_cols = [c for c in physchem_cols if c in hits_df.columns]

    # Natural sort order of conditions
    try:
        condition_order = sorted(
            hits_df["condition"].unique(),
            key=lambda x: int(x.split('_')[0][1:])
        )
    except Exception as e:
        logger.warning("[plot_del_hits] Failed to parse numeric condition order: %s", e)
        condition_order = sorted(hits_df["condition"].unique())

    hits_df["condition"] = pd.Categorical(
        hits_df["condition"],
        categories=condition_order,
        ordered=True
    )

    # Heatmap
    df_unique = hits_df.groupby(["NsynthonID", "condition"], as_index=False)[plot_col].max()
    top_synthons = (
        df_unique.groupby("NsynthonID")[plot_col].max().sort_values(ascending=False).head(heatmap_top_n).index
    )
    df_top = df_unique[df_unique["NsynthonID"].isin(top_synthons)]
    heat_df = df_top.pivot(index="NsynthonID", columns="condition", values=plot_col).fillna(0)
    heat_df = heat_df.loc[heat_df.max(axis=1).sort_values(ascending=False).index]

    fig_height = max(6, len(heat_df) * 0.25)
    plt.figure(figsize=(12, fig_height))
    sns.heatmap(heat_df, cmap="berlin", linewidths=0.5)
    plt.title(f"{plot_col} per synthon across conditions (Top {heatmap_top_n})")
    plt.ylabel("NsynthonID")
    plt.xlabel("Condition")
    plt.tight_layout()
    heatmap_file = output_dir / "heatmap_p_active.png"
    plt.savefig(heatmap_file, dpi=300)
    plt.close()
    plot_files["heatmap"] = heatmap_file

    # Stripplot 
    plt.figure(figsize=(12, 6))
    sns.stripplot(data=hits_df, x="condition", y=plot_col, jitter=True, size=5, alpha=0.7, order=condition_order)
    plt.ylabel(plot_col)
    plt.title(f"Distribution of {plot_col} per condition")
    plt.xlabel("Condition")
    plt.tight_layout()
    strip_file = output_dir / "stripplot_p_active.png"
    plt.savefig(strip_file, dpi=300)
    plt.close()
    plot_files["stripplot"] = strip_file

    # Top-x barplot across conditions
    top_synthons = hits_df.groupby("NsynthonID")[plot_col].max().sort_values(ascending=False).head(top_n).index
    plot_df = hits_df[hits_df["NsynthonID"].isin(top_synthons)].copy()
    synthon_order = plot_df.groupby("NsynthonID")[plot_col].max().sort_values(ascending=False).index

    plt.figure(figsize=(max(10, 1.2 * len(synthon_order)), 6))
    sns.barplot(
        data=plot_df,
        x="NsynthonID",
        y=plot_col,
        hue="condition",
        order=synthon_order,
        hue_order=condition_order
    )
    plt.xticks(rotation=90)
    plt.ylabel(plot_col)
    plt.xlabel("NsynthonID")
    plt.title(f"Top {top_n} synthons across conditions by {plot_col}")
    plt.legend(title="Condition")
    plt.tight_layout()
    bar_file = output_dir / "top_hits_barplot.png"
    plt.savefig(bar_file, dpi=300)
    plt.close()
    plot_files["barplot"] = bar_file

    logger.debug(
        "[plot_del_hits] hits_df rows: %d | full_hits_df rows: %d | unique synthons (full): %d",
        len(hits_df),
        len(full_hits_df),
        full_hits_df["NsynthonID"].nunique()
    )

    # UpSet plot 
    # Use the ordered list of conditions for columns
    conds = full_hits_df["condition"].unique()

    # Natural sort
    condition_order = sorted(conds, key=lambda x: int(x.split('_')[0][1:]))

    indicator_df = pd.DataFrame(
        False,
        index=full_hits_df["NsynthonID"].unique(),
        columns=condition_order
    )

    for cond in condition_order:
        synthons_in_cond = full_hits_df.loc[full_hits_df["condition"] == cond, "NsynthonID"].unique()
        indicator_df.loc[synthons_in_cond, cond] = True

    upset_data = from_indicators(indicator_df.columns, indicator_df)
    plt.figure(figsize=(10, 6))
    UpSet(upset_data, show_counts=True, sort_by="cardinality").plot()
    plt.suptitle("Overlap of hits across conditions")
    upset_file = output_dir / "hits_upset.png"
    plt.savefig(upset_file, dpi=300)
    plt.close()
    plot_files["upset"] = upset_file

    # Ridgeline plots vs physchem
    valid_physchem = [c for c in physchem_cols if c in full_hits_df.columns and full_hits_df[c].notna().sum() >= 2]

    if not valid_physchem:
        logger.warning("[plot_del_hits] Ridgeline plots skipped: no valid numeric physchem columns with >=2 values found.")
    else:
        logger.info("[plot_del_hits] Generating ridgeline plots for columns: %s", valid_physchem)

    for col in valid_physchem:
        ridge_df = full_hits_df.copy()
        ridge_df = ridge_df[ridge_df[col].notna()]

        condition_order = sorted(
            ridge_df["condition"].unique(),
            key=lambda x: int(x.split('_')[0][1:]),
            # default=[]
        )

        plt.figure(figsize=(10, max(6, 0.5*len(condition_order))))

        for i, cond in enumerate(condition_order):
            vals = ridge_df.loc[ridge_df["condition"] == cond, col].values
            if len(vals) < 2 or np.all(vals == vals[0]):
                logger.debug(f"[plot_del_hits] Skipping condition {cond} for {col}: insufficient or constant values")
                continue

            # Use Gaussian KDE
            kde = gaussian_kde(vals, bw_method=0.5)
            x = np.linspace(vals.min(), vals.max(), 200)
            y = kde(x)

            # Offset y by row index for ridgeline
            plt.fill_between(x, y + i, i, alpha=0.7)
            plt.plot(x, y + i, color="black", lw=0.5)
            plt.text(x[0], i + 0.05, cond, fontsize=8)

        plt.xlabel(col)
        plt.ylabel("Condition")
        plt.yticks([])  # remove all y-axis ticks
        plt.title(f"Ridgeline plot of {col} per condition")
        plt.tight_layout()
        ridge_file = Path(output_dir) / f"ridgeline_{col}.png"
        plt.savefig(ridge_file, dpi=300)
        plt.close()
        plot_files[f"ridgeline_{col}"] = ridge_file

    # Scatterplots vs physchem
    valid_physchem = [c for c in physchem_cols if c in hits_df.columns and hits_df[c].notna().any()]

    if not valid_physchem:
        logger.warning("[plot_del_hits] Scatterplots skipped: no valid numeric physchem columns with non-NaN values found.")
    else:
        logger.info("[plot_del_hits] Scatterplots will be generated for columns: %s", valid_physchem)
        
        for col in valid_physchem:
            jitter_x = (hits_df[col].max() - hits_df[col].min()) * 0.01
            jitter_y = (hits_df[plot_col].max() - hits_df[plot_col].min()) * 0.01
            x = hits_df[col] + np.random.normal(0, jitter_x, size=len(hits_df))
            y = hits_df[plot_col] + np.random.normal(0, jitter_y, size=len(hits_df))

            g = sns.FacetGrid(hits_df.assign(x_jitter=x, y_jitter=y),
                            col="condition", col_wrap=3, height=4, sharex=False, sharey=False)
            g.map_dataframe(sns.scatterplot, x="x_jitter", y="y_jitter", alpha=0.4, s=20)
            g.set_axis_labels(col, plot_col)
            g.set_titles(col_template="{col_name}")
            g.add_legend()

            plt.tight_layout()
            scatter_file = output_dir / f"{plot_col}_vs_{col}_facet.png"
            plt.savefig(scatter_file, dpi=300)
            plt.close()
            plot_files[f"scatter_{col}"] = scatter_file


    return plot_files

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




