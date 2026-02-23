import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm
from tqdm_joblib import tqdm_joblib
from matplotlib.colors import ListedColormap
from collections import Counter
from math import log2

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

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

def get_cdr_regions_from_mutable(mutable_positions):
    """
    Converts a sorted list of mutable positions into contiguous regions.
    Returns dict with keys CDR1, CDR2, ... and values (start, end)
    """
    if not mutable_positions:
        return {}

    mutable_positions = sorted(mutable_positions)
    regions = []
    start = mutable_positions[0]
    prev = mutable_positions[0]

    for pos in mutable_positions[1:]:
        if pos != prev + 1:
            regions.append((start, prev))
            start = pos
        prev = pos
    regions.append((start, prev))

    # Give names automatically: CDR1, CDR2, ...
    cdr_regions = {f"CDR{i+1}": (s, e) for i, (s, e) in enumerate(regions)}
    return cdr_regions

def generate_summary_plots(
    df_all,
    centroid_history,
    ancestor_colour_map,
    context,
    data_dir,
    plots_dir,
    ancestor_to_seq,
    aa_to_int,
    cdr_regions,
    framework_positions=None,
    max_variants=150,
    n_jobs=10
):
    """
    Generate all summary plots and per-ancestor visualisations.
    Parallelizes per-ancestor plots for speed.
    """

    logger.info("Generating summary plots and ancestor prioritisation CSV...")

    # ------------------------------------------------------------------
    # 1. Centroid trajectories
    # ------------------------------------------------------------------
    plt.figure(figsize=(8, 7))
    for ancestor_id, history in centroid_history.items():
        if len(history) < 2:
            continue
        hx = [h["umap1"] for h in history]
        hy = [h["umap2"] for h in history]
        rounds = [h["round"] for h in history]
        plt.plot(hx, hy, marker="o", linewidth=2, alpha=0.9,
                 color=ancestor_colour_map.get(ancestor_id, "gray"), label=f"Ancestor {ancestor_id}")
        plt.scatter(hx[-1], hy[-1], s=140, edgecolor="black", zorder=5,
                    color=ancestor_colour_map.get(ancestor_id, "gray"))
        plt.text(hx[-1], hy[-1], f"R{rounds[-1]}", fontsize=9)

    plt.title("Pareto Centroid Trajectories Across Rounds")
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.legend(title="Ancestor ID", bbox_to_anchor=(1.05, 1))
    plt.tight_layout()
    plt.savefig(plots_dir / "ancestor_centroid_trajectories.png", dpi=150)
    plt.close()

    # ------------------------------------------------------------------
    # 2. Ancestor prioritisation CSV
    # ------------------------------------------------------------------
    ancestor_summary = (
        df_all
        .groupby("ancestor_id")
        .agg(
            total_descendants=("sequence", "count"),
            pareto_rate=("pareto_flag", "mean"),
            best_score=("score", "max"),
            mean_score=("score", "mean"),
            max_mutations=("mut_count", "max"),
            survived_rounds=("round", "max"),
        )
        .reset_index()
        .sort_values("best_score", ascending=False)
    )
    ancestor_summary.to_csv(data_dir / "ancestor_prioritisation_summary.csv", index=False)
    logger.info("Ancestor prioritisation summary saved")

    # ------------------------------------------------------------------
    # 3. Round-over-round performance (physchem + CI + Pareto)
    # ------------------------------------------------------------------

    df_stats = pd.DataFrame(context["round_stats"])

    props = [
        c.replace("median_", "")
        for c in df_stats.columns
        if c.startswith("median_") and not c.endswith("_lci") and not c.endswith("_uci")
    ]

    n_props = len(props)

    fig, axes = plt.subplots(
        n_props, 1,
        figsize=(8, 4 * n_props),
        sharex=True
    )

    if n_props == 1:
        axes = [axes]

    for ax, prop in zip(axes, props):

        med_col = f"median_{prop}"
        lci_col = f"{med_col}_lci"
        uci_col = f"{med_col}_uci"
        pareto_col = f"pareto_median_{prop}"

        # Main median line
        ax.plot(
            df_stats["round"],
            df_stats[med_col],
            marker="o",
            label="All candidates"
        )

        # CI shading
        if lci_col in df_stats and uci_col in df_stats:
            ax.fill_between(
                df_stats["round"],
                df_stats[lci_col],
                df_stats[uci_col],
                alpha=0.25
            )

        # Pareto overlay
        if pareto_col in df_stats:
            ax.plot(
                df_stats["round"],
                df_stats[pareto_col],
                linestyle="--",
                marker="s",
                label="Pareto only"
            )

        ax.set_ylabel(prop)
        ax.set_title(f"{prop} across rounds")
        ax.grid(True)
        ax.legend()

    axes[-1].set_xlabel("Round")

    plt.tight_layout()
    plt.savefig(plots_dir / "round_over_round_all_properties.png", dpi=150)
    plt.close()

    # ------------------------------------------------------------------
    # 4. Per-ancestor plot families
    # ------------------------------------------------------------------
    residue_dir   = plots_dir / "ancestor_residue_heatmaps"
    chemistry_dir = plots_dir / "ancestor_chemistry_heatmaps"
    entropy_dir   = plots_dir / "ancestor_entropy_tracks"

    for d in [residue_dir, chemistry_dir, entropy_dir]:
        d.mkdir(parents=True, exist_ok=True)

    ancestors = sorted(df_all["ancestor_id"].unique())
    logger.info(f"Generating per-ancestor plots for {len(ancestors)} ancestors using {n_jobs if n_jobs != -1 else 'all'} cores")

    def plot_ancestor(anc):
    # Residue identity heatmap 
        plot_residue_annotated_heatmap(
            df_all=df_all,
            ancestor_to_seq=ancestor_to_seq,
            ancestor_id=anc,
            aa_to_int=aa_to_int,
            cdr_regions=cdr_regions,
            framework_positions=framework_positions,
            max_variants=max_variants,
            save_path=residue_dir / f"ancestor_{anc}_residue_heatmap.png"
        )

        # Chemistry heatmap
        plot_chemistry_annotated_heatmap(
            df_all=df_all,
            ancestor_to_seq=ancestor_to_seq,
            ancestor_id=anc,
            cdr_regions=cdr_regions,
            framework_positions=framework_positions,
            max_variants=max_variants,
            save_path=chemistry_dir / f"ancestor_{anc}_chemistry_heatmap.png"
        )

        # Mutation entropy track
        plot_mutation_entropy_track_for_ancestor(
            df_all=df_all,
            ancestor_id=anc,
            ancestor_to_seq=ancestor_to_seq,
            cdr_regions=cdr_regions,
            framework_positions=framework_positions,
            max_variants=max_variants,
            save_dir=entropy_dir
        )

    residue_dir   = plots_dir / "ancestor_residue_heatmaps"
    chemistry_dir = plots_dir / "ancestor_chemistry_heatmaps"
    entropy_dir   = plots_dir / "ancestor_entropy_tracks"
    for d in [residue_dir, chemistry_dir, entropy_dir]:
        d.mkdir(parents=True, exist_ok=True)

    ancestors = sorted(df_all["ancestor_id"].unique())

    # Plotting (parallel)
    with tqdm_joblib(tqdm(desc="Ancestor plots", total=len(ancestors))) as progress_bar:
        Parallel(n_jobs=-1)(
            delayed(plot_ancestor)(anc)
            for anc in ancestors
        )

    logger.info("All summary plots and ancestor visualisations complete")

def mutation_matrix_against_parent(parent_seq, seqs):
    parent = np.array(list(parent_seq))
    mat = np.array([parent != np.array(list(seq)) for seq in seqs], dtype=int)
    return mat

def residue_matrix_against_parent(parent_seq, seqs, aa_to_int):
    parent = np.array(list(parent_seq))
    mat = []
    for seq in seqs:
        row = [-1 if aa == parent[p] else aa_to_int[aa] for p, aa in enumerate(seq)]
        mat.append(row)
    return np.array(mat)

def plot_residue_annotated_heatmap(
    df_all,
    ancestor_to_seq,
    ancestor_id,
    aa_to_int,
    cdr_regions,
    framework_positions=None,
    max_variants=150,
    save_path=None
):
    """
    Heatmap showing mutations against parent sequence for a single ancestor.
    """
    df_a = (
        df_all[df_all["ancestor_id"] == ancestor_id]
        .sort_values("score", ascending=False)
        .head(max_variants)
    )
    if df_a.empty:
        print(f"No sequences found for ancestor {ancestor_id}")
        return

    parent = ancestor_to_seq[ancestor_id]
    mut_mask = mutation_matrix_against_parent(parent, df_a["sequence"].tolist())
    res_mat = residue_matrix_against_parent(parent, df_a["sequence"], aa_to_int)

    fig, ax = plt.subplots(figsize=(18, max(4, len(df_a) * 0.05)))
    fig.subplots_adjust(top=0.85)

    # Background: mutation mask
    ax.imshow(mut_mask, aspect="auto", cmap="Greys", vmin=0, vmax=1, interpolation="nearest")

    # Foreground: residue identity
    masked_res = np.ma.masked_where(res_mat < 0, res_mat)
    aa_list = list(aa_to_int.keys())
    aa_cmap = ListedColormap(sns.color_palette("tab20", len(aa_list)))
    aa_cmap.set_bad(color=(0, 0, 0, 0))  # transparent
    im = ax.imshow(masked_res, aspect="auto", cmap=aa_cmap,
                   vmin=-0.5, vmax=len(aa_list) - 0.5, interpolation="nearest")

    # Framework shading
    if framework_positions:
        for start, end in framework_positions:
            ax.axvspan(start - 0.5, end - 0.5, color="lightgray", alpha=0.25, zorder=3)

    # CDR annotations
    for label, (start, end) in cdr_regions.items():
        ax.axvline(start - 0.5, color="black", linewidth=1)
        ax.axvline(end - 0.5, color="black", linewidth=1)
        mid = (start + end) / 2
        ax.text(mid, 1.05, label, ha="center", va="bottom",
                fontsize=11, fontweight="bold", transform=ax.get_xaxis_transform())

    ax.set_xlabel("Residue position")
    ax.set_ylabel("Variants")
    ax.set_title(f"Ancestor {ancestor_id}: mutation + residue identity")

    # Colorbar for residues
    cbar = plt.colorbar(im, ax=ax, ticks=range(len(aa_list)))
    cbar.ax.set_yticklabels(aa_list)
    cbar.set_label("Mutated residue")

    # Make room for labels above
    ax.set_ylim(ax.get_ylim()[0] + 2.5, ax.get_ylim()[1])
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150)
        plt.close()
    else:
        plt.show()

AA_CHEMISTRY = {
    # Hydrophobic
    "A": "hydrophobic", "V": "hydrophobic", "I": "hydrophobic",
    "L": "hydrophobic", "M": "hydrophobic", "F": "hydrophobic",
    "W": "hydrophobic", "Y": "hydrophobic",

    # Polar uncharged
    "S": "polar (S, T, N, Q)", "T": "polar (S, T, N, Q)", "N": "polar (S, T, N, Q)", "Q": "polar (S, T, N, Q)",

    # Positive
    "K": "positive (K, R, H)", "R": "positive (K, R, H)", "H": "positive (K, R, H)",

    # Negative
    "D": "negative (D, E)", "E": "negative (D, E)",

    # Special
    "C": "Other (C, G, P)", "G": "Other (C, G, P)", "P": "Other (C, G, P)",
}

chem_groups = ["hydrophobic", "polar (S, T, N, Q)", "positive (K, R, H)", "negative (D, E)", "Other (C, G, P)"]
chem_to_int = {c: i for i, c in enumerate(chem_groups)}

def chemistry_matrix_against_parent(parent_seq, seqs):
    """
    Returns matrix (N, L):
    -1 = same as parent
     0..len(chem_groups)-1 = chemistry class
    """
    parent = np.array(list(parent_seq))
    mat = []

    for seq in seqs:
        row = []
        for p, aa in enumerate(seq):
            if aa == parent[p]:
                row.append(-1)
            else:
                row.append(chem_to_int[AA_CHEMISTRY[aa]])
        mat.append(row)

    return np.array(mat)

chem_colors = {
    "hydrophobic": "#FDB462", 
    "polar (S, T, N, Q)": "#80B1D3",
    "positive (K, R, H)": "#FB8072",
    "negative (D, E)": "#B3DE69",
    "Other (C, G, P)": "#BC80BD",
}

chem_cmap = ListedColormap(
    [chem_colors[c] for c in chem_groups]
)
chem_cmap.set_bad(color=(0, 0, 0, 0))  # transparent for no mutation

def plot_chemistry_annotated_heatmap(
    *,
    df_all,
    ancestor_id,
    ancestor_to_seq,
    cdr_regions,
    framework_positions=None,
    max_variants=150,
    save_path=None,
):
    """
    Chemistry-annotated mutation heatmap for a single ancestor.
    """

    df_a = (
        df_all[df_all["ancestor_id"] == ancestor_id]
        .sort_values("score", ascending=False)
        .head(max_variants)
    )

    if df_a.empty:
        logger.warning(f"No sequences found for ancestor {ancestor_id}")
        return

    parent = ancestor_to_seq[ancestor_id]
    mut_mask = mutation_matrix_against_parent(parent, df_a["sequence"].tolist())
    chem_mat = chemistry_matrix_against_parent(parent, df_a["sequence"].tolist())

    fig, ax = plt.subplots(figsize=(18, max(4, len(df_a) * 0.05)))

    # Mutation background 
    ax.imshow(
        mut_mask,
        aspect="auto",
        cmap="Greys",
        vmin=0,
        vmax=1,
        interpolation="nearest"
    )

    # Chemistry overlay
    masked_chem = np.ma.masked_where(chem_mat < 0, chem_mat)
    im = ax.imshow(
        masked_chem,
        aspect="auto",
        cmap=chem_cmap,
        vmin=-0.5,
        vmax=len(chem_groups) - 0.5,
        interpolation="nearest"
    )

    # Framework shading
    if framework_positions:
        for start, end in framework_positions:
            ax.axvspan(
                start - 0.5,
                end - 0.5,
                color="lightgray",
                alpha=0.15,
                zorder=3
            )

    # CDR annotations
    for label, (start, end) in cdr_regions.items():
        ax.axvline(start - 0.5, color="black", linewidth=1)
        ax.axvline(end - 0.5, color="black", linewidth=1)
        mid = (start + end) / 2
        ax.text(
            mid,
            1.05,
            label,
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold",
            transform=ax.get_xaxis_transform()
        )

    ax.set_xlabel("Residue position")
    ax.set_ylabel("Variants")
    ax.set_title(f"Ancestor {ancestor_id}: mutation chemistry map")

    # Chemistry legend
    cbar = plt.colorbar(im, ax=ax, ticks=range(len(chem_groups)))
    cbar.ax.set_yticklabels(chem_groups)
    cbar.set_label("Mutated residue chemistry")

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150)
        plt.close()
    else:
        plt.show()

def mutation_entropy_against_parent(parent_seq, seqs):
    """
    Returns array (L,) of mutation entropy per position.
    Entropy computed ONLY over mutated residues.
    """
    parent = list(parent_seq)
    L = len(parent)
    entropy = np.zeros(L)

    for pos in range(L):
        muts = [seq[pos] for seq in seqs if seq[pos] != parent[pos]]

        if len(muts) == 0:
            entropy[pos] = 0.0
            continue

        counts = Counter(muts)
        total = sum(counts.values())
        H = 0.0
        for c in counts.values():
            p = c / total
            H -= p * log2(p)
        entropy[pos] = H

    return entropy


def plot_mutation_entropy_track_for_ancestor(
    df_all,
    ancestor_id,
    ancestor_to_seq,
    cdr_regions=None,
    framework_positions=None,
    max_variants=150,
    save_dir="ancestor_entropy_tracks"
):
    """
    Compute mutation entropy track for a single ancestor and save the plot.
    """
    df_a = (
        df_all[df_all["ancestor_id"] == ancestor_id]
        .sort_values("score", ascending=False)
        .head(max_variants)
    )

    if df_a.empty:
        print(f"No sequences found for ancestor {ancestor_id}")
        return

    parent_seq = ancestor_to_seq[ancestor_id]
    seqs = df_a["sequence"].tolist()

    entropy = mutation_entropy_against_parent(parent_seq, seqs)

    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(16, 4))
    plt.bar(range(len(entropy)), entropy, color="steelblue")

    if framework_positions:
        for start, end in framework_positions:
            plt.axvspan(
                start - 0.5,
                end - 0.5,
                color="lightgray",
                alpha=0.15
            )

    if cdr_regions:
        for name, (start, end) in cdr_regions.items():
            plt.axvspan(start, end, alpha=0.15, color="orange")
            plt.text((start + end) / 2, max(entropy) * 1.05,
                     name, ha="center", fontsize=10, fontweight="bold")

    plt.ylabel("Entropy (bits)")
    plt.xlabel("Residue position")
    plt.title(f"Ancestor {ancestor_id}: mutation entropy")
    plt.tight_layout()

    save_path = save_dir / f"ancestor_{ancestor_id}_entropy.png"
    plt.savefig(save_path, dpi=150)
    plt.close()
    return save_path

def plot_learning_curves_per_property(df_metrics, out_path):
    """
    Creates executive summary:
    Rows = properties
    Cols = [Train Spearman | Train RMSE | Future Spearman]
    """
    properties = sorted(df_metrics["property"].unique())
    n_props = len(properties)

    fig, axes = plt.subplots(
        n_props, 3,
        figsize=(15, 4 * n_props),
        sharex=True
    )

    if n_props == 1:
        axes = axes.reshape(1, -1)

    for i, prop in enumerate(properties):
        dfp = df_metrics[df_metrics["property"] == prop]

        # Train Spearman
        ax = axes[i, 0]
        ax.plot(dfp["round"], dfp["train_spearman"], marker="o")
        ax.fill_between(
            dfp["round"],
            dfp["train_spearman_low"],
            dfp["train_spearman_high"],
            alpha=0.25
        )
        ax.set_ylim(-1, 1)
        ax.set_title(f"{prop}\nTrain Spearman")
        ax.axhline(0, color="gray", linestyle="--")

        # Train RMSE
        ax = axes[i, 1]
        ax.plot(dfp["round"], dfp["train_rmse"], marker="o")
        ax.fill_between(
            dfp["round"],
            dfp["train_rmse_low"],
            dfp["train_rmse_high"],
            alpha=0.25
        )
        ax.set_title(f"{prop}\nTrain RMSE")

        # Future Spearman
        ax = axes[i, 2]
        ax.plot(dfp["round"], dfp["future_spearman"], marker="o")
        ax.fill_between(
            dfp["round"],
            dfp["future_spearman_low"],
            dfp["future_spearman_high"],
            alpha=0.25
        )
        ax.set_ylim(-1, 1)
        ax.set_title(f"{prop}\nFuture-only Spearman")
        ax.axhline(0, color="gray", linestyle="--")

    for ax in axes[-1, :]:
        ax.set_xlabel("Active learning round")

    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()