import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

def plot_del_hits(hits_csv: str, output_dir: str = "outputs/plots", top_n: int = 10):
    """
    Generate heatmap, strip plot, and bar plot from DEL hits CSV.

    Args:
        hits_csv: Path to DEL hits CSV (full per-condition data)
        output_dir: Directory to save plots
        top_n: Number of top hits per condition for bar plot
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(hits_csv)

    # --- Heatmap ---
    heat_df = df.pivot(index="NsynthonID", columns="condition", values="p_active")
    plt.figure(figsize=(12, 8))
    sns.heatmap(heat_df, cmap="viridis", linewidths=0.5)
    plt.title("Per-synthon p_active across conditions")
    plt.ylabel("NsynthonID")
    plt.xlabel("Condition")
    plt.tight_layout()
    heatmap_file = output_dir / "heatmap_p_active.png"
    plt.savefig(heatmap_file, dpi=300)
    plt.close()
    print(f"Heatmap saved to {heatmap_file}")

    # --- Strip Plot ---
    plt.figure(figsize=(12, 6))
    sns.stripplot(data=df, x="condition", y="p_active", jitter=True, size=5, alpha=0.7)
    plt.title("Distribution of p_active per condition")
    plt.ylabel("p_active")
    plt.xlabel("Condition")
    plt.tight_layout()
    strip_file = output_dir / "stripplot_p_active.png"
    plt.savefig(strip_file, dpi=300)
    plt.close()
    print(f"Strip plot saved to {strip_file}")

    # --- Top Hits Bar Plot ---
    top_df = (
        df.sort_values("p_active", ascending=False)
        .groupby("condition")
        .head(top_n)
        .copy()
    )
    plt.figure(figsize=(14, 6))
    sns.barplot(data=top_df, x="NsynthonID", y="p_active", hue="condition")
    plt.xticks(rotation=90)
    plt.title(f"Top {top_n} synthons per condition")
    plt.tight_layout()
    bar_file = output_dir / "top_hits_barplot.png"
    plt.savefig(bar_file, dpi=300)
    plt.close()
    print(f"Top hits bar plot saved to {bar_file}")
