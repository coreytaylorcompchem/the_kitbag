import torch
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from pathlib import Path
from torch_geometric.loader import DataLoader

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def evaluate(context, config):
    """
    Evaluation for CYP3A4 GIN model.
    Computes predictions, metrics, and saves a true vs predicted scatter plot
    with RMSE, MAE, and R² included in the plot title.
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

    y_true = []
    y_pred = []

    with torch.no_grad():
        for batch in val_loader:
            batch = batch.to(device)

            if hasattr(batch, 'global_features'):
                global_feat = batch.global_features.to(device).float()
                if global_feat.dim() == 1:
                    global_feat = global_feat.unsqueeze(0).repeat(batch.num_graphs, 1)
                batch.global_features = global_feat

            batch.y = torch.nan_to_num(batch.y, nan=0.0).float()

            pred = model(batch)
            y_true.extend(batch.y.cpu().numpy())
            y_pred.extend(pred.cpu().numpy())

    # Compute metrics
    mse = mean_squared_error(y_true, y_pred)
    rmse = np.sqrt(mse)
    mae = mean_absolute_error(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)

    logger.info(f"\nRMSE: {rmse:.4f}")
    logger.info(f"MAE: {mae:.4f}")
    logger.info(f"R²: {r2:.4f}")

    # Save plot
    plot_dir = Path(config.get("eval_plot_dir", "_images"))
    plot_dir.mkdir(parents=True, exist_ok=True)
    plot_path = plot_dir / "true_vs_predicted_pic50s.png"

    plt.figure(figsize=(12, 6))
    plt.scatter(y_true, y_pred, alpha=0.5, label="Predictions")
    plt.plot([min(y_true), max(y_true)], [min(y_true), max(y_true)], 'r--', label="Ideal")
    plt.xlabel("True pIC50")
    plt.ylabel("Predicted pIC50")
    plt.title(f"CYP3A4 pIC50 Prediction (GIN Model)\nRMSE={rmse:.4f}, MAE={mae:.4f}, R²={r2:.4f}")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300)
    if config.get("save_plots", True):
        plt.show()
    else:
        plt.close()

    return {
        "evaluation_metrics": {"MSE": mse, "RMSE": rmse, "MAE": mae, "R2": r2},
        "eval_plot_path": str(plot_path)
    }
