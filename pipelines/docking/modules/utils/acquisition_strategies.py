import numpy as np
from sklearn.metrics import pairwise_distances


def select_acquisition_candidates(preds, uncertainty, batch_size, strategy, ligand_indices, ligand_fps=None):
    strategy = strategy.lower()

    if strategy == "uncertainty_sampling":# Select top ligands with highest prediction uncertainty (exploitation + exploration).
        selected = np.argsort(-uncertainty)[:batch_size]

    elif strategy == "greedy": # Pick ligands with lowest predicted docking scores (pure exploitation).
        selected = np.argsort(preds)[:batch_size]

    elif strategy == "random": # Select randomly (baseline).
        selected = np.random.choice(len(preds), size=batch_size, replace=False)

    elif strategy == "expected_improvement": # Select based on improvement over current best, using mean & std.
        best_score = np.min(preds)  # Lower is better
        improvement = best_score - preds
        ei = improvement * uncertainty
        selected = np.argsort(-ei)[:batch_size]

    elif strategy == "thompson_sampling": # Select based on improvement over current best, using mean & std.
        selected = np.random.choice(len(preds), size=batch_size, replace=False, p=softmax(-preds))

    elif strategy == "diverse_uncertainty": # Add MaxMin fingerprint diversity on top of uncertainty sampling.
        if ligand_fps is None:
            raise ValueError("Ligand fingerprints are required for diverse acquisition.")
        selected = diverse_uncertainty_selection(preds, uncertainty, ligand_fps, batch_size)

    else:
        raise ValueError(f"Unknown acquisition strategy: {strategy}")

    return [ligand_indices[i] for i in selected]


def softmax(x):
    e_x = np.exp(x - np.max(x))
    return e_x / e_x.sum()


def diverse_uncertainty_selection(preds, uncertainty, fps, batch_size):
    """
    MaxMin diversity selection using uncertainty as primary criterion.
    """
    scores = uncertainty
    sorted_indices = np.argsort(-scores)
    selected = [sorted_indices[0]]

    while len(selected) < batch_size:
        remaining = list(set(sorted_indices) - set(selected))
        dists = pairwise_distances(fps[remaining], fps[selected], metric='jaccard')
        min_dists = np.min(dists, axis=1)
        next_idx = remaining[np.argmax(min_dists)]
        selected.append(next_idx)

    return selected
