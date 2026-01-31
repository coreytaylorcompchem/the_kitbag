import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

def plot_label_histograms(
    df_before,
    df_after,
    label_col="pIC50",
    out_dir="_images",
    bins=40,
    title_suffix="",
):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(14, 5))

    # Before
    plt.subplot(1, 2, 1)
    plt.hist(df_before[label_col], bins=bins, alpha=0.8)
    plt.title(f"Before sampling{title_suffix}")
    plt.xlabel(label_col)
    plt.ylabel("Count")
    plt.grid(True)

    # After
    plt.subplot(1, 2, 2)
    plt.hist(df_after[label_col], bins=bins, alpha=0.8)
    plt.title(f"After sampling{title_suffix}")
    plt.xlabel(label_col)
    plt.ylabel("Count")
    plt.grid(True)

    plt.tight_layout()

    out_path = out_dir / "label_distribution_before_after.png"
    plt.savefig(out_path, dpi=300)
    plt.close()

    return str(out_path)

# compute distance by seed
def hamming_distance(s1, s2):
    return sum(1 for a,b in zip(s1,s2) if a != b)

# generate one-hot embeddings
def one_hot_encode(seq, seq_len, aa_list, aa_to_idx):
    oh = np.zeros((seq_len, len(aa_list)))
    for i, aa in enumerate(seq):
        oh[i, aa_to_idx[aa]] = 1
    return oh.flatten()

def compute_multi_condition_pareto(df, properties):
    """
    Identify non-dominated candidates across all specified properties
    Returns boolean array indicating Pareto-optimal rows
    """
    data = df[properties].values
    is_pareto = np.ones(data.shape[0], dtype=bool)
    for i, row in enumerate(data):
        if is_pareto[i]:
            is_dominated = np.all(data >= row, axis=1) & np.any(data > row, axis=1)
            is_pareto[is_dominated] = False
    return is_pareto

def plot_round_radar(df_candidates, round_idx, save_path, properties):
    """
    Create a radar plot comparing Pareto vs non-Pareto sequences for multi-condition properties
    """
    # Normalise props 0-1
    df_norm = df_candidates.copy()
    for prop in properties:
        min_val = df_norm[prop].min()
        max_val = df_norm[prop].max()
        df_norm[prop] = (df_norm[prop] - min_val) / (max_val - min_val)

    # Compute mean (pareto and non-pareto)
    mean_pareto = df_norm[df_norm["pareto_flag"]][properties].mean()
    mean_non_pareto = df_norm[~df_norm["pareto_flag"]][properties].mean()

    labels = properties
    N = len(labels)
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
    angles += angles[:1]
    
    fig, ax = plt.subplots(figsize=(6,6), subplot_kw=dict(polar=True))

    # Pareto mean
    values = mean_pareto.tolist()
    values += values[:1]
    ax.plot(angles, values, color='gold', linewidth=2, label='Pareto')
    ax.fill(angles, values, color='gold', alpha=0.25)

    # Non-Pareto mean
    values = mean_non_pareto.tolist()
    values += values[:1]
    ax.plot(angles, values, color='gray', linewidth=2, label='Non-Pareto')
    ax.fill(angles, values, color='gray', alpha=0.25)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_yticks([0.0, 0.5, 1.0])
    ax.set_ylim(0,1)
    ax.set_title(f"Round {round_idx} Property Profile (Radar)", pad=20)
    ax.legend(loc='upper right', bbox_to_anchor=(1.2,1.1))
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    plt.close()