import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

def plot_performance(per_cycle_df, output_dir, sample_size=None, batch_size=None, baseline_df=None):
    plt.figure(figsize=(10, 6))

    # Active learning plot
    sns.boxplot(x="cycle", y="score", data=per_cycle_df, palette="Blues", showfliers=False, width=0.5)
    sns.stripplot(x="cycle", y="score", data=per_cycle_df, color='black', alpha=0.4, jitter=True, size=3)

    # Baseline plot
    if baseline_df is not None and not baseline_df.empty:
        sns.boxplot(x="cycle", y="score", data=baseline_df, palette="Reds", showfliers=False, width=0.5)
        sns.stripplot(x="cycle", y="score", data=baseline_df, color='darkred', alpha=0.3, jitter=True, size=3)

    # Title with optional parameters
    title = "Docking Score Distribution Over Active Learning Cycles"
    if sample_size is not None or batch_size is not None:
        title += f"\n(Sample size = {sample_size}, Batch size = {batch_size})"
    plt.title(title)

    plt.ylabel("Docking Score (Lower is Better)")
    plt.xlabel("Cycle")
    plt.grid(True)
    plt.tight_layout()

    # Create custom legend
    legend_elements = [
        Patch(facecolor='lightblue', edgecolor='black', label='Active Learning'),
        Patch(facecolor='lightcoral', edgecolor='black', label='Random Baseline'),
        Line2D([0], [0], marker='o', color='w', label='Individual Docking Scores',
               markerfacecolor='black', markersize=5, alpha=0.5)
    ]
    plt.legend(handles=legend_elements, loc='upper right')

    # Save plot
    output_path = Path(output_dir) / "docking_score_by_cycle.png"
    plt.savefig(output_path)
    plt.close()
