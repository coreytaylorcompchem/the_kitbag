import numpy as np
import torch

def denormalise(y: np.ndarray, mean: np.ndarray, std: np.ndarray) -> np.ndarray:
    """
    Denormalise an array of predictions or targets from normalized space back to original scale.
    
    Args:
        y: np.ndarray of shape (N, outputs) - normalized data
        mean: np.ndarray of shape (outputs,) - mean used for normalization
        std: np.ndarray of shape (outputs,) - std dev used for normalization

    Returns:
        np.ndarray - denormalised data
    """
    return y * std + mean

def normalise_graph_labels(graphs):
    """Normalize graph.y labels feature-wise across dataset."""
    all_y = np.stack([g.y.squeeze(0).numpy() for g in graphs])
    y_mean = all_y.mean(axis=0)
    y_std = all_y.std(axis=0) + 1e-8

    for g in graphs:
        g.y = torch.tensor(
            ((g.y.squeeze(0).numpy() - y_mean) / y_std),
            dtype=torch.float32
        ).unsqueeze(0)

    return y_mean, y_std
