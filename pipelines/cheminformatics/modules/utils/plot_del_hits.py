import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

def plot_del_hits(hits_df: pd.DataFrame, output_dir: str = "outputs/plots", top_n: int = 10):
    """
    Generate three plots from DEL hit dataframe:
    1. Heatmap of p_active per synthon x condition
    2. Strip plot of p_active per condition
    3. Bar plot of top hits per condition
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Heatmap
    heat_df = hits_df.pivot(index="NsynthonID", columns="condition", values="p_active")
    plt.figure(figsize=(12, 8))
    sns.heatmap(heat_df, cmap="viridis", linewidths=0.5)
    plt.title("Per-synthon p_active across conditions")
    plt.ylabel("NsynthonID")
    plt.xlabel("Condition")
    plt.tight_layout()
    heatmap_file = output_dir / "heatmap_p_active.png"
    plt.savefig(heatmap_file, dpi=300)
    plt.close()

    # Strip plot
    plt.figure(figsize=(12, 6))
    sns.stripplot(data=hits_df, x="condition", y="p_active", jitter=True, size=5, alpha=0.7)
    plt.title("Distribution of p_active per condition")
    plt.ylabel("p_active")
    plt.xlabel("Condition")
    plt.tight_layout()
    strip_file = output_dir / "stripplot_p_active.png"
    plt.savefig(strip_file, dpi=300)
    plt.close()

    # Top hits bar plot
    top_df = hits_df.sort_values("p_active", ascending=False).groupby("condition").head(top_n)
    plt.figure(figsize=(14, 6))
    sns.barplot(data=top_df, x="NsynthonID", y="p_active", hue="condition")
    plt.xticks(rotation=90)
    plt.title(f"Top {top_n} synthons per condition")
    plt.tight_layout()
    bar_file = output_dir / "top_hits_barplot.png"
    plt.savefig(bar_file, dpi=300)
    plt.close()

    return {
        "heatmap": heatmap_file,
        "stripplot": strip_file,
        "barplot": bar_file
    }
