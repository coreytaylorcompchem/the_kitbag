import pandas as pd
import numpy as np
import seaborn as sns

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec

from pathlib import Path
from upsetplot import from_indicators, UpSet
from scipy.stats import gaussian_kde
from itertools import combinations

from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

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
    physchem_filters: dict | None = None,  # NEW
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

    reduced_df = hits_df            # reduced (UMAP-sampled)
    full_df = full_hits_df.copy()   # complete hits
    full_df["condition"] = full_df["condition"].astype(str)

    # ----------------------------
    # Unique synthon dataset
    # ----------------------------
    if "NsynthonID" in full_df.columns:
        syn_df = full_df.drop_duplicates("NsynthonID").copy()
    else:
        syn_df = full_df.copy()

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

    # # Heatmap
    # df_unique = reduced_df.groupby(["NsynthonID", "condition"], as_index=False)[plot_col].max()
    # top_synthons = (
    #     df_unique.groupby("NsynthonID")[plot_col].max().sort_values(ascending=False).head(heatmap_top_n).index
    # )
    # df_top = df_unique[df_unique["NsynthonID"].isin(top_synthons)]
    # heat_df = df_top.pivot(index="NsynthonID", columns="condition", values=plot_col).fillna(0)
    # heat_df = heat_df.loc[heat_df.max(axis=1).sort_values(ascending=False).index]

    # fig_height = max(6, len(heat_df) * 0.25)
    # plt.figure(figsize=(12, fig_height))
    # sns.heatmap(heat_df, cmap="vanimo", linewidths=0.5)
    # plt.title(f"{plot_col} per synthon across conditions (Top {heatmap_top_n})")
    # plt.ylabel("NsynthonID")
    # plt.xlabel("Condition")
    # plt.tight_layout()
    # heatmap_file = output_dir / "heatmap_p_active.png"
    # plt.savefig(heatmap_file, dpi=300)
    # plt.close()
    # plot_files["heatmap"] = heatmap_file

    # Heatmap (clustered)

    df_unique = reduced_df.groupby(
        ["NsynthonID", "condition"],
        as_index=False
    )[plot_col].max()

    top_synthons = (
        df_unique.groupby("NsynthonID")[plot_col]
        .max()
        .sort_values(ascending=False)
        .head(heatmap_top_n)
        .index
    )

    df_top = df_unique[df_unique["NsynthonID"].isin(top_synthons)]

    heat_df = df_top.pivot(
        index="NsynthonID",
        columns="condition",
        values=plot_col
    ).fillna(0)

    # -----------------------------
    # Row-normalise (pattern clustering)
    # -----------------------------
    heat_norm = heat_df.sub(
        heat_df.mean(axis=1), axis=0
    ).div(
        heat_df.std(axis=1).replace(0, 1), axis=0)

    # -----------------------------
    # Cluster only if possible
    # -----------------------------
    if heat_norm.shape[0] > 1:

        row_dist = pdist(heat_norm.values, metric="euclidean")
        row_linkage = linkage(row_dist, method="ward")
        row_order = leaves_list(row_linkage)

        heat_df_clustered = heat_df.iloc[row_order]

    else:
        # Not enough rows to cluster
        heat_df_clustered = heat_df.copy()

    # -----------------------------
    # Plot
    # -----------------------------
    fig_height = max(6, len(heat_df_clustered) * 0.25)

    plt.figure(figsize=(12, fig_height))
    sns.heatmap(
        heat_df_clustered,
        cmap="vanimo",
        linewidths=0.5
    )

    plt.title(f"{plot_col} per synthon across conditions (Clustered Top {heatmap_top_n})")
    plt.ylabel("NsynthonID")
    plt.xlabel("Condition")
    plt.tight_layout()

    heatmap_file = output_dir / "heatmap_clustered.png"
    plt.savefig(heatmap_file, dpi=300)
    plt.close()

    plot_files["heatmap"] = heatmap_file

    # ----------------------------
    # Box + strip plot by condition
    # ----------------------------

    plt.figure(figsize=(12, 6))

    # Boxplot for summary stats
    sns.boxplot(
        data=full_df,
        x="condition",
        y=plot_col,
        order=condition_order,
        whis=1.5,         
        color="lightgray",
        fliersize=0    
    )

    # Overlay stripplot for individual points
    sns.stripplot(
        data=full_df,
        x="condition",
        y=plot_col,
        size=5,
        jitter=True,
        alpha=0.2,
        order=condition_order,
        color="dodgerblue"
    )

    plt.ylabel(plot_col)
    plt.xlabel("Condition")
    plt.title(f"Distribution of {plot_col} per condition.")
    plt.tight_layout()

    box_strip_file = output_dir / "box_stripplot_p_active.png"
    plt.savefig(box_strip_file, dpi=300)
    plt.close()

    plot_files["box_stripplot"] = box_strip_file

    logger.info("[plot_del_hits] Box + strip plot saved: %s", box_strip_file)

    # ----------------------------
    # UMAP colored by SMARTS severity
    # ----------------------------
    if {"umap_x", "umap_y", "smarts_severity"}.issubset(syn_df.columns):

        plt.figure(figsize=(7,6))

        sns.scatterplot(
            data=syn_df,
            x="umap_x",
            y="umap_y",
            hue="smarts_severity",
            palette={
                "GREEN": "seagreen",
                "ORANGE": "orange",
                "RED": "crimson"
            },
            s=30,
            alpha=0.85
        )

        plt.title("Chemical space colored by SMARTS severity")
        plt.xlabel("UMAP 1")
        plt.ylabel("UMAP 2")

        plt.legend(title="SMARTS severity", frameon=False)

        plt.tight_layout()

        sev_umap_file = output_dir / "umap_smarts_severity.png"
        plt.savefig(sev_umap_file, dpi=300)
        plt.close()

        plot_files["umap_smarts_severity"] = sev_umap_file

    else:
        logger.warning("[plot_del_hits] Missing columns for SMARTS UMAP plot")

    # ----------------------------
    # SMARTS composition per cluster
    # ----------------------------
    if {"cluster_id", "smarts_severity"}.issubset(syn_df.columns):

        cluster_sev = (
            syn_df
            .groupby(["cluster_id", "smarts_severity"])
            .size()
            .unstack(fill_value=0)
        )

        cluster_sev = cluster_sev.div(cluster_sev.sum(axis=1), axis=0)

        cluster_sev = cluster_sev.reindex(columns=["GREEN","ORANGE","RED"], fill_value=0)

        cluster_sev.plot(
            kind="bar",
            stacked=True,
            figsize=(8,6),
            color={
                "GREEN": "seagreen",
                "ORANGE": "orange",
                "RED": "crimson"
            }
        )

        plt.ylabel("Fraction of synthons")
        plt.xlabel("UMAP cluster")
        plt.title("SMARTS composition per chemical cluster")

        plt.legend(title="SMARTS severity", frameon=False)

        plt.tight_layout()

        cluster_smarts_file = output_dir / "cluster_smarts_composition.png"
        plt.savefig(cluster_smarts_file, dpi=300)
        plt.close()

        plot_files["cluster_smarts_composition"] = cluster_smarts_file

    else:
        logger.warning("[plot_del_hits] Missing columns for cluster SMARTS plot")


    # Top-x barplot across conditions
    top_synthons = (
        reduced_df.groupby("NsynthonID")[plot_col]
        .max()
        .sort_values(ascending=False)
        .head(top_n)
        .index
    )

    plot_df = reduced_df[reduced_df["NsynthonID"].isin(top_synthons)].copy()
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
    if full_hits_df.empty:
        logger.warning("[plot_del_hits] UpSet plot skipped: no hits available.")
    else:

        conds = full_hits_df["condition"].unique()

        try:
            condition_order = sorted(conds, key=lambda x: int(x.split('_')[0][1:]))
        except Exception:
            condition_order = sorted(conds)

        indicator_df = pd.DataFrame(
            False,
            index=full_hits_df["NsynthonID"].unique(),
            columns=condition_order
        )

        for cond in condition_order:
            synthons_in_cond = full_hits_df.loc[
                full_hits_df["condition"] == cond,
                "NsynthonID"
            ].unique()

            indicator_df.loc[synthons_in_cond, cond] = True

        if indicator_df.empty:
            logger.warning("[plot_del_hits] UpSet plot skipped: empty indicator matrix.")
        else:
            upset_data = from_indicators(indicator_df.columns, indicator_df)

            plt.figure(figsize=(10, 6))
            UpSet(upset_data, show_counts=True, sort_by="cardinality").plot()
            plt.suptitle("Overlap of hits across conditions")

            upset_file = output_dir / "hits_upset.png"
            plt.savefig(upset_file, dpi=300)
            plt.close()

            plot_files["upset"] = upset_file

    # ----------------------------
    # Ridgeline plots vs physchem
    # ----------------------------
    valid_physchem = [c for c in physchem_cols if c in full_hits_df.columns and full_hits_df[c].notna().sum() >= 2]

    if not valid_physchem:
        logger.warning("[plot_del_hits] Ridgeline plots skipped: no valid numeric physchem columns with >=2 values found.")
    else:
        logger.info("[plot_del_hits] Generating ridgeline plots for columns: %s", valid_physchem)

    for col in valid_physchem:
        ridge_df = full_hits_df.copy()

        # Apply physchem filters ONLY for ridgeline plots
        if physchem_filters:
            for fcol, (low, high) in physchem_filters.items():
                if fcol in ridge_df.columns:
                    before = len(ridge_df)
                    ridge_df = ridge_df[ridge_df[fcol].between(low, high)]
                    logger.debug(
                        "[plot_del_hits][ridgeline] Applied %s in [%.2f, %.2f]: %d → %d rows",
                        fcol, low, high, before, len(ridge_df)
                    )

        ridge_df = ridge_df[ridge_df[col].notna()]

        # Only keep conditions with enough points
        condition_counts = ridge_df.groupby("condition")[col].count()
        condition_order = sorted(
            condition_counts[condition_counts >= 5].index,
            key=lambda x: int(x.split('_')[0][1:])
        )

        if not condition_order:
            logger.warning(f"[plot_del_hits] Ridgeline plot skipped for {col}: no conditions with >=5 hits")
            continue

        plt.figure(figsize=(10, max(6, 0.5*len(condition_order))))

        # x-axis range based on actual values
        x_min, x_max = ridge_df[col].min(), ridge_df[col].max()
        x_range = np.linspace(x_min, x_max, 200)

        for i, cond in enumerate(condition_order):
            vals = ridge_df.loc[ridge_df["condition"] == cond, col].values
            if len(vals) < 2 or np.all(vals == vals[0]):
                logger.debug(f"[plot_del_hits] Skipping condition {cond} for {col}: insufficient or constant values")
                continue

            # Gaussian KDE with adaptive bandwidth (otherwise Ridge plots can squish to the left or right)
            bw = max(0.1, min(0.5, 1.0 / np.sqrt(len(vals))))
            kde = gaussian_kde(vals, bw_method=bw)
            y = kde(x_range)

            # Normalize per condition
            y /= y.max()

            # Offset y by row index
            plt.fill_between(x_range, i, y + i, alpha=0.7)
            plt.plot(x_range, y + i, color="black", lw=0.5)
            plt.text(x_min, i + 0.05, cond, fontsize=8)

        plt.xlabel(col)
        plt.ylabel("Condition")
        plt.yticks([])  # remove y-axis ticks
        plt.title(f"Ridgeline plot of {col} per condition")
        plt.tight_layout()

        ridge_file = Path(output_dir) / f"ridgeline_{col}.png"
        plt.savefig(ridge_file, dpi=300)
        plt.close()
        plot_files[f"ridgeline_{col}"] = ridge_file


    # Scatterplots vs physchem
    valid_physchem = [c for c in physchem_cols if c in full_df.columns and full_df[c].notna().any()]

    if not valid_physchem:
        logger.warning("[plot_del_hits] Scatterplots skipped: no valid numeric physchem columns with non-NaN values found.")
    else:
        logger.info("[plot_del_hits] Scatterplots will be generated for columns: %s", valid_physchem)
        
        for col in valid_physchem:
            jitter_x = (full_df[col].max() - full_df[col].min()) * 0.01
            jitter_y = (full_df[plot_col].max() - full_df[plot_col].min()) * 0.01

            x = full_df[col] + np.random.normal(0, jitter_x, size=len(full_df))
            y = full_df[plot_col] + np.random.normal(0, jitter_y, size=len(full_df))

            g = sns.FacetGrid(
                full_df.assign(x_jitter=x, y_jitter=y),
                col="condition",
                col_wrap=3,
                height=4,
                sharex=False,
                sharey=False
            )
            g.map_dataframe(sns.scatterplot, x="x_jitter", y="y_jitter", alpha=0.4, s=20)
            g.set_axis_labels(col, plot_col)
            g.set_titles(col_template="{col_name}")
            g.add_legend()

            plt.tight_layout()
            scatter_file = output_dir / f"{plot_col}_vs_{col}_facet.png"
            plt.savefig(scatter_file, dpi=300)
            plt.close()
            plot_files[f"scatter_{col}"] = scatter_file
    
    # ----------------------------
    # Faceted UMAP by condition
    # ----------------------------
    if {"umap_x", "umap_y"}.issubset(reduced_df.columns):

        conditions = sorted(reduced_df["condition"].unique())
        n_cols = 3
        n_rows = int(np.ceil(len(conditions) / n_cols))

        fig = plt.figure(figsize=(4 * n_cols, 4 * n_rows + 1))
        gs = GridSpec(n_rows + 1, n_cols, height_ratios=[4]*n_rows + [0.2], hspace=0.3)

        vmin = reduced_df[plot_col].min()
        vmax = reduced_df[plot_col].max()

        axes = []
        for i, cond in enumerate(conditions):
            row = i // n_cols
            col = i % n_cols
            ax = fig.add_subplot(gs[row, col])
            axes.append(ax)

            # background points
            ax.scatter(
                reduced_df["umap_x"],
                reduced_df["umap_y"],
                s=10,
                c="lightgray",
                alpha=0.2,
                linewidth=0
            )

            # highlight points for this condition
            sub = reduced_df[reduced_df["condition"] == cond]
            sc = ax.scatter(
                sub["umap_x"],
                sub["umap_y"],
                c=sub[plot_col],
                cmap="plasma_r",
                vmin=vmin,
                vmax=vmax,
                s=20,
                alpha=0.9,
                linewidth=0
            )

            ax.set_title(cond)
            ax.set_xticks([])
            ax.set_yticks([])

        # turn off unused axes
        for i in range(len(conditions), n_rows*n_cols):
            ax = fig.add_subplot(gs[i // n_cols, i % n_cols])
            ax.axis("off")

        # colorbar in the bottom row of GridSpec
        cbar_ax = fig.add_subplot(gs[-1, :])
        cbar = fig.colorbar(sc, cax=cbar_ax, orientation="horizontal", label=plot_col)
        cbar_ax.xaxis.set_ticks_position('bottom')

        fig.suptitle("UMAP of DEL hits (condition-highlighted)", y=1.02)
        umap_file = output_dir / "umap_by_condition.png"
        plt.savefig(umap_file, dpi=300, bbox_inches="tight")
        plt.close()
        plot_files["umap"] = umap_file

    else:
        logger.warning(
            "[plot_del_hits] UMAP scatter skipped: umap_x / umap_y not found in dataframe"
        )
    
    # ----------------------------
    # Library representation: Full vs UMAP-reduced
    # ----------------------------
    if {"NsynthonID"}.issubset(full_df.columns):

        def extract_library(series):
            return series.str.split("-", n=1).str[0]

        full_unique = full_df.drop_duplicates("NsynthonID")
        reduced_unique = reduced_df.drop_duplicates("NsynthonID")

        full_lib = extract_library(full_unique["NsynthonID"])
        reduced_lib = extract_library(reduced_unique["NsynthonID"])

        # Proportions
        full_counts = full_lib.value_counts(normalize=True)
        reduced_counts = reduced_lib.value_counts(normalize=True)

        # Sort by full representation
        full_counts = full_counts.sort_values(ascending=False)

        # Keep top N for readability
        top_n_libs = 15
        top_libs = full_counts.head(top_n_libs).index

        full_top = full_counts.reindex(top_libs).fillna(0)
        reduced_top = reduced_counts.reindex(top_libs).fillna(0)

        # Add "Other"
        full_other = full_counts.drop(top_libs).sum()
        reduced_other = reduced_counts.drop(top_libs).sum()

        full_top.loc["Other"] = full_other
        reduced_top.loc["Other"] = reduced_other

        lib_plot_df = pd.DataFrame({
            "library": full_top.index,
            "Full hits": full_top.values,
            "Reduced (UMAP)": reduced_top.values
        })

        lib_plot_df = lib_plot_df.sort_values("Full hits", ascending=True)

        plt.figure(figsize=(10, 8))

        bar_width = 0.4
        y_pos = np.arange(len(lib_plot_df))

        plt.barh(
            y_pos - bar_width/2,
            lib_plot_df["Full hits"],
            height=bar_width,
            label="Full hits",
            alpha=0.8
        )

        plt.barh(
            y_pos + bar_width/2,
            lib_plot_df["Reduced (UMAP)"],
            height=bar_width,
            label="Reduced (UMAP)",
            alpha=0.8
        )

        plt.yticks(y_pos, lib_plot_df["library"])
        plt.xlabel("Proportion of hits")
        plt.title("Library representation before and after UMAP reduction")
        plt.legend()
        plt.tight_layout()

        lib_file = output_dir / "library_full_vs_reduced.png"
        plt.savefig(lib_file, dpi=300)
        plt.close()

        plot_files["library_bias"] = lib_file
    
    return plot_files


def plot_del_qc_metrics(
    df: pd.DataFrame,
    output_dir: str = "outputs/plots_qc",
):
    """
    Generate DEL QC plots to analyse quality of evidence:
    - Effect size vs uncertainty
    - HDI overlap with zero
    - Pre vs post count enrichment
    - Hit consistency across conditions
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    qc_files = {}
    df = df.copy()

    # ----------------------------
    # Derived QC quantities
    # ----------------------------

    df["hdi_width"] = df["beta_hdi_upper"] - df["beta_hdi_lower"]
    df["hdi_crosses_zero"] = (
        (df["beta_hdi_lower"] <= 0) & (df["beta_hdi_upper"] >= 0)
    )

    # ----------------------------
    # Plot 1: Effect size vs uncertainty
    # ----------------------------

    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        data=df,
        x="beta_mean",
        y="hdi_width",
        hue="condition",
        alpha=0.6,
        s=40
    )

    plt.axvline(0, color="black", lw=1, linestyle="--")
    plt.xlabel("Effect size (beta mean)")
    plt.ylabel("Uncertainty (HDI width)")
    plt.title("Effect size vs uncertainty per synthon")
    plt.tight_layout()

    f = output_dir / "qc_effect_vs_uncertainty.png"
    plt.savefig(f, dpi=300)
    plt.close()
    qc_files["effect_vs_uncertainty"] = f

    # ----------------------------
    # Plot 2: HDI overlap with zero
    # ----------------------------

    hdi_summary = (
        df.assign(
            hdi_status=np.where(
                df["beta_hdi_upper"] < 0, "HDI < 0",
                np.where(df["beta_hdi_lower"] > 0, "HDI > 0", "HDI overlaps 0")
            )
        )
        .groupby("hdi_status")
        .size()
        .reset_index(name="count")
    )

    plt.figure(figsize=(6, 4))
    sns.barplot(
        data=hdi_summary,
        x="hdi_status",
        y="count",
        order=["HDI > 0", "HDI overlaps 0", "HDI < 0"],
        palette="muted"
    )
    plt.ylabel("Number of synthons")
    plt.xlabel("")
    plt.title("HDI overlap with zero (screen-level signal)")
    plt.tight_layout()

    f = output_dir / "qc_hdi_overlap.png"
    plt.savefig(f, dpi=300)
    plt.close()
    qc_files["hdi_overlap"] = f

    # ----------------------------
    # Plot 3: Posterior enrichment vs effect size
    # ----------------------------

    if {"beta_mean", "enrichment_mean"}.issubset(df.columns):
        plot_df = df.copy()

        if not plot_df.empty:
            plt.figure(figsize=(7, 6))
            sc = plt.scatter(
                np.log10(plot_df["enrichment_mean"] + 1e-6),  # avoid log(0)
                plot_df["beta_mean"],
                c=plot_df["p_active"],
                cmap="coolwarm",
                alpha=0.6,
                s=40
            )

            plt.axhline(0, color="black", lw=1, linestyle="--")
            plt.xlabel("log10(posterior enrichment mean)")
            plt.ylabel("Effect size (beta_mean)")
            plt.title("Posterior enrichment vs effect size.")

            cbar = plt.colorbar(sc)
            cbar.set_label("p_active")

            plt.tight_layout()
            f = output_dir / "qc_posterior_enrichment_vs_effect.png"
            plt.savefig(f, dpi=300)
            plt.close()

            qc_files["posterior_enrichment_vs_effect"] = f
    else:
        logger.info(
            "[QC] posterior enrichment not present — skipping enrichment vs effect plot"
        )

    # ----------------------------
    # Plot 4: Hit consistency across conditions
    # ----------------------------

    active_threshold = 0.5

    hit_counts = (
        df[df["p_active"] >= active_threshold]
        .groupby("NsynthonID")["condition"]
        .nunique()
        .reset_index(name="n_conditions")
    )

    plt.figure(figsize=(6, 4))
    sns.histplot(
        hit_counts["n_conditions"],
        bins=range(1, hit_counts["n_conditions"].max() + 2),
        discrete=True
    )

    plt.xlabel("Number of conditions with activity")
    plt.ylabel("Number of synthons")
    plt.title(f"Hit consistency across conditions (p_active ≥ {active_threshold})")
    plt.tight_layout()

    f = output_dir / "qc_hit_consistency.png"
    plt.savefig(f, dpi=300)
    plt.close()
    qc_files["hit_consistency"] = f

    return qc_files
