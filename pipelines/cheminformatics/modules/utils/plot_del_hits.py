import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from upsetplot import from_indicators, UpSet
import warnings

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

    # --- Ensure correct types ---
    hits_df = hits_df.copy()
    hits_df["p_active"] = pd.to_numeric(hits_df["p_active"], errors="coerce")
    hits_df["condition"] = hits_df["condition"].astype(str)

    # --- Identify numeric physchem columns ---
    if physchem_cols is None:
        physchem_cols = [
            c for c in hits_df.columns
            if c not in {"NsynthonID", "condition", "beta_mean", "beta_hdi_lower", "beta_hdi_upper", "p_active"}
            and pd.api.types.is_numeric_dtype(hits_df[c])
        ]

    for col in physchem_cols:
        hits_df[col] = pd.to_numeric(hits_df[col], errors="coerce")

    # --- Heatmap (Top synthons by max p_active) ---
    # Take max p_active per synthon × condition to remove duplicates
    df_unique = hits_df.groupby(["NsynthonID", "condition"], as_index=False)["p_active"].max()

    # Select top N synthons by max p_active across all conditions
    top_synthons = df_unique.groupby("NsynthonID")["p_active"].max() \
                            .sort_values(ascending=False).head(heatmap_top_n).index
    df_top = df_unique[df_unique["NsynthonID"].isin(top_synthons)]

    # Pivot for heatmap
    heat_df = df_top.pivot(index="NsynthonID", columns="condition", values="p_active").fillna(0)

    # Sort synthons by max p_active (so top synthon at top)
    heat_df = heat_df.loc[heat_df.max(axis=1).sort_values(ascending=False).index]

    # Adjust figure height dynamically
    fig_height = max(6, len(heat_df) * 0.25)  # 0.25 inch per synthon
    plt.figure(figsize=(12, fig_height))
    sns.heatmap(heat_df, cmap="viridis", linewidths=0.5)
    plt.title(f"Per-synthon p_active across conditions (Top {heatmap_top_n})")
    plt.ylabel("NsynthonID")
    plt.xlabel("Condition")
    plt.tight_layout()
    heatmap_file = output_dir / "heatmap_p_active.png"
    plt.savefig(heatmap_file, dpi=300)
    plt.close()
    plot_files["heatmap"] = heatmap_file

    # --- Stripplot ---
    plt.figure(figsize=(12, 6))
    sns.stripplot(data=hits_df, x="condition", y="p_active", jitter=True, size=5, alpha=0.7)
    plt.title("Distribution of p_active per condition")
    plt.ylabel("p_active")
    plt.xlabel("Condition")
    plt.tight_layout()
    strip_file = output_dir / "stripplot_p_active.png"
    plt.savefig(strip_file, dpi=300)
    plt.close()
    plot_files["stripplot"] = strip_file

    # --- Top hits bar plot ---
    top_df = hits_df.sort_values("p_active", ascending=False).groupby("condition").head(top_n)
    plt.figure(figsize=(14, 6))
    sns.barplot(data=top_df, x="NsynthonID", y="p_active", hue="condition")
    plt.xticks(rotation=90)
    plt.title(f"Top {top_n} synthons per condition")
    plt.tight_layout()
    bar_file = output_dir / "top_hits_barplot.png"
    plt.savefig(bar_file, dpi=300)
    plt.close()
    plot_files["barplot"] = bar_file

    # --- UpSet plot ---
    # Create indicator matrix for hits
    conds = hits_df["condition"].unique()
    indicator_df = pd.DataFrame(False, index=hits_df["NsynthonID"].unique(), columns=conds)
    for cond in conds:
        indicator_df.loc[hits_df.loc[hits_df["condition"] == cond, "NsynthonID"], cond] = True

    upset_data = from_indicators(indicator_df.columns, indicator_df)
    plt.figure(figsize=(10, 6))
    UpSet(upset_data, show_counts=True, sort_by="cardinality").plot()
    plt.suptitle("Overlap of hits across conditions")
    upset_file = output_dir / "hits_upset.png"
    plt.savefig(upset_file, dpi=300)
    plt.close()
    plot_files["upset"] = upset_file

    # --- Scatterplots vs physchem ---
    physchem_cols = [c for c in physchem_cols if c in hits_df.columns and pd.api.types.is_numeric_dtype(hits_df[c])]
    if physchem_cols:
        n_cols = 3
        n_rows = (len(physchem_cols) + n_cols - 1) // n_cols
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*5, n_rows*4))
        if n_rows * n_cols == 1:
            axes = [axes]
        else:
            axes = axes.flatten()

        for i, col in enumerate(physchem_cols):
            sns.scatterplot(
                data=hits_df,
                x=col,
                y="p_active",
                hue="condition",
                ax=axes[i],
                alpha=0.7,
                palette="tab10"
            )
            axes[i].set_title(f"p_active vs {col}")

        # Remove empty subplots
        for j in range(i+1, len(axes)):
            fig.delaxes(axes[j])

        # Single legend
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper right', title="Condition")

        plt.tight_layout(rect=[0,0,0.85,1])
        scatter_file = output_dir / "p_active_vs_physchem.png"
        plt.savefig(scatter_file, dpi=300)
        plt.close()
        plot_files["scatter"] = scatter_file

    return plot_files
