import torch
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error, mean_absolute_error
from scipy.stats import pearsonr
from pathlib import Path
from torch_geometric.loader import DataLoader

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def evaluate(context, config):
    """
    Evaluation for CYP3A4 GIN model.
    Computes predictions, metrics, and saves a true vs predicted scatter plot
    with RMSE, MAE, and Pearson r included in the plot title.
    """

    model = context.get("trained_model") or context.get("model")
    if model is None:
        raise ValueError("No trained model found in context.")

    device = context.get("device", torch.device("cpu"))
    model.to(device)
    model.eval()

    val_loader = context.get("val_loader")
    if val_loader is None:
        raise ValueError("No validation DataLoader found in context.")

    # Collect predictions and true labels
    y_true, y_pred = [], []

    with torch.no_grad():
        for batch in val_loader:
            batch = batch.to(device)

            if hasattr(batch, "global_features"):
                global_feat = batch.global_features.to(device).float()
                if global_feat.dim() == 1:
                    global_feat = global_feat.unsqueeze(0).repeat(batch.num_graphs, 1)
                batch.global_features = global_feat

            pred = model(batch)

            # Flatten both arrays to 1D
            y_true.append(batch.y.view(-1).cpu())
            y_pred.append(pred.view(-1).cpu())

    # Concatenate to single 1D arrays
    y_true = torch.cat(y_true).numpy()
    y_pred = torch.cat(y_pred).numpy()

    # Mask out invalid values
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    y_true = y_true[mask]
    y_pred = y_pred[mask]

    # Compute metrics
    mse = mean_squared_error(y_true, y_pred)
    rmse = np.sqrt(mse)
    mae = mean_absolute_error(y_true, y_pred)
    r, _ = pearsonr(y_true, y_pred)  # Pearson correlation

    logger.info(f"\nRMSE: {rmse:.4f}")
    logger.info(f"MAE: {mae:.4f}")
    logger.info(f"Pearson r: {r:.4f}")

    # Save scatter plot
    plot_dir = Path(config.get("eval_plot_dir", "_images"))
    plot_dir.mkdir(parents=True, exist_ok=True)
    plot_path = plot_dir / "true_vs_predicted_pic50s.png"

    plt.figure(figsize=(12, 6))
    plt.scatter(y_true, y_pred, alpha=0.5, label="Predictions")
    plt.plot([min(y_true), max(y_true)], [min(y_true), max(y_true)], 'r--', label="Ideal")
    plt.xlabel("True pIC50")
    plt.ylabel("Predicted pIC50")
    plt.title(f"CYP3A4 pIC50 Prediction (GIN Model)\nRMSE={rmse:.4f}, MAE={mae:.4f}, Pearson r={r:.4f}")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300)
    if config.get("save_plots", True):
        plt.show()
    else:
        plt.close()

    return {
        "evaluation_metrics": {"MSE": mse, "RMSE": rmse, "MAE": mae, "Pearson_r": r},
        "eval_plot_path": str(plot_path)
    }
