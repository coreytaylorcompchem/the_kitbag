import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from upsetplot import from_indicators, UpSet
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def plot_del_hits(hits_df: pd.DataFrame, output_dir: str = "outputs/plots",
                  top_n: int = 10, physchem_cols=None, heatmap_top_n: int = 100):
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

    # Choose plotting column
    plot_col = "p_active_cond" if "p_active_cond" in hits_df.columns else "p_active"
    hits_df[plot_col] = pd.to_numeric(hits_df[plot_col], errors="coerce")

    # --- Determine numeric physchem columns ---
    if physchem_cols is None:
        physchem_cols = [
            c for c in hits_df.columns
            if c not in {"NsynthonID", "condition", "beta_mean", "beta_hdi_lower", "beta_hdi_upper", "p_active", "p_active_cond"}
            and pd.api.types.is_numeric_dtype(hits_df[c])
        ]
    physchem_cols = [c for c in physchem_cols if c in hits_df.columns]

    # --- Natural order of conditions ---
    try:
        condition_order = sorted(
            hits_df["condition"].unique(),
            key=lambda x: int(x.split('_')[0][1:])  # 'C2_E3' -> 2
        )
    except Exception as e:
        logger.warning("[plot_del_hits] Failed to parse numeric condition order: %s", e)
        condition_order = sorted(hits_df["condition"].unique())

    hits_df["condition"] = pd.Categorical(
        hits_df["condition"],
        categories=condition_order,
        ordered=True
    )

    # --- Heatmap ---
    df_unique = hits_df.groupby(["NsynthonID", "condition"], as_index=False)[plot_col].max()
    top_synthons = (
        df_unique.groupby("NsynthonID")[plot_col].max().sort_values(ascending=False).head(heatmap_top_n).index
    )
    df_top = df_unique[df_unique["NsynthonID"].isin(top_synthons)]
    heat_df = df_top.pivot(index="NsynthonID", columns="condition", values=plot_col).fillna(0)
    heat_df = heat_df.loc[heat_df.max(axis=1).sort_values(ascending=False).index]

    fig_height = max(6, len(heat_df) * 0.25)
    plt.figure(figsize=(12, fig_height))
    sns.heatmap(heat_df, cmap="viridis", linewidths=0.5)
    plt.title(f"{plot_col} per synthon across conditions (Top {heatmap_top_n})")
    plt.ylabel("NsynthonID")
    plt.xlabel("Condition")
    plt.tight_layout()
    heatmap_file = output_dir / "heatmap_p_active.png"
    plt.savefig(heatmap_file, dpi=300)
    plt.close()
    plot_files["heatmap"] = heatmap_file

    # --- Stripplot ---
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

    # --- Top-X barplot across conditions ---
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

    # --- UpSet plot ---
    # Use the ordered list of conditions for columns
    conds = hits_df["condition"].unique()
    # Preserve your natural order, e.g., C1_E3, C2_E3, ...
    condition_order = sorted(conds, key=lambda x: int(x.split('_')[0][1:]))

    # Create indicator matrix with ordered columns
    indicator_df = pd.DataFrame(False, index=hits_df["NsynthonID"].unique(), columns=condition_order)
    for cond in condition_order:
        indicator_df.loc[hits_df.loc[hits_df["condition"] == cond, "NsynthonID"], cond] = True

    # Build UpSet plot (combinations still sorted by size)
    upset_data = from_indicators(indicator_df.columns, indicator_df)
    plt.figure(figsize=(10, 6))
    UpSet(upset_data, show_counts=True, sort_by="cardinality").plot()
    plt.suptitle("Overlap of hits across conditions")
    upset_file = output_dir / "hits_upset.png"
    plt.savefig(upset_file, dpi=300)
    plt.close()
    plot_files["upset"] = upset_file

    # --- Scatterplots vs physchem ---
    valid_physchem = [c for c in physchem_cols if hits_df[c].notna().any()]
    if not valid_physchem:
        logger.warning("[plot_del_hits] Scatterplots skipped: no valid numeric physchem columns with non-NaN values found.")
    else:
        n_cols = 3
        n_rows = (len(valid_physchem) + n_cols - 1) // n_cols
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*5, n_rows*4), squeeze=False)
        axes = axes.flatten()

        for i, col in enumerate(valid_physchem):
            sns.scatterplot(data=hits_df, x=col, y=plot_col, hue="condition", ax=axes[i], alpha=0.7)
            axes[i].set_title(f"{plot_col} vs {col}")

        for j in range(i+1, len(axes)):
            fig.delaxes(axes[j])

        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, title="Condition", loc="upper right")
        plt.tight_layout(rect=[0, 0, 0.85, 1])

        scatter_file = output_dir / f"{plot_col}_vs_physchem.png"
        plt.savefig(scatter_file, dpi=300)
        plt.close()
        plot_files["scatter"] = scatter_file

    return plot_files
