import os
import itertools
import matplotlib.pyplot as plt
import pickle
import numpy as np
from pathlib import Path

import torch
from torch_geometric.loader import DataLoader
from torch.utils.data import random_split

from sklearn.metrics import mean_squared_error, r2_score
from sklearn.linear_model import LinearRegression

from models.pose_ranker import PoseScoringGNN

from pipeline.task_registry import register_task

from modules.utils.graph_construction import load_all_targets

from models.utils.data_utils import normalise_graph_labels
from models.utils.data_utils import denormalise
from models.utils.train_utils import train_epoch, eval_epoch, EarlyStopping

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)


@register_task("load_targets_prepare_graphs", category="Pose ranker", description="Load MD trajs and parameterise as graphs.")
def load_targets_prepare_graphs(config: dict, context: dict):
    """
    Load MD trajectories, featurise ligand + context, and produce per-frame graphs.
    Config:
        input_trajs: str - base directory with target subfolders
        ligand_resname: str - residue name of ligand (default: UNK)
        save_graphs_to: str - output pickle file
    """
    base_dir = Path(config.get("input_trajs", ""))
    ligand_resname = config.get("ligand_resname", "UNK")
    out_path = Path(config.get("save_graphs_to", "graphs.pkl"))

    if not base_dir.exists():
        raise FileNotFoundError(f"Input trajectory directory not found: {base_dir}")

    logger.info(f"Loading and featurising trajectories from {base_dir} (ligand={ligand_resname})")
    graphs = load_all_targets(str(base_dir), ligand_resname=ligand_resname)

    logger.info(f"Writing {len(graphs)} graphs to {out_path}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "wb") as f:
        pickle.dump(graphs, f)

    context["graphs_path"] = str(out_path)
    return {"graphs": graphs, "graphs_path": str(out_path)}

@register_task("load_model_spec", category="Pose ranker", description="Load Pose Ranker model spec.")
def load_model_spec(config: dict, context: dict):
    """
    Instantiate the PoseScoringGNN model and add to context.

    Config:
        in_dim: int - number of input features (required)
        hidden_dim: int - hidden layer size (default 128)
        n_outputs: int - output size (default 9)
        device: str - 'cpu' or 'cuda' (default: 'cpu')
    """
    in_dim = config.get("in_dim")
    if in_dim is None:
        # Try to infer from graphs if available
        graphs = context.get("graphs")
        if graphs and len(graphs) > 0:
            in_dim = graphs[0].x.shape[1]
        else:
            raise ValueError("in_dim not provided and cannot infer from graphs")

    hidden_dim = config.get("hidden_dim", 128)
    n_outputs = config.get("n_outputs", 9)
    device = torch.device(config.get("device", "cpu"))

    model = PoseScoringGNN(in_dim, hidden_dim=hidden_dim, n_outputs=n_outputs).to(device)
    context["model"] = model
    context["device"] = device

    logger.info(f"Model initialised: PoseScoringGNN(in_dim={in_dim}, hidden_dim={hidden_dim}, n_outputs={n_outputs}) on {device}")

    return {"model": model, "device": device}

@register_task("prepare_dataloaders", category="Pose ranker", description="Normalise labels and create train/val DataLoaders")
def prepare_dataloaders(config: dict, context: dict):
    graphs = context.get("graphs")
    if not graphs:
        raise ValueError("No graphs found in context. Run data loading task first.")

    batch_size = config.get("batch_size", 32)
    val_split = config.get("val_split", 0.2)
    shuffle = config.get("shuffle", True)
    normalize_labels_flag = config.get("normalize_labels", True)

    if normalize_labels_flag:
        y_mean, y_std = normalise_graph_labels(graphs)
    else:
        y_mean = np.zeros(graphs[0].y.shape[-1])
        y_std = np.ones(graphs[0].y.shape[-1])

    n_total = len(graphs)
    val_size = int(val_split * n_total)
    train_size = n_total - val_size

    train_dataset, val_dataset = random_split(graphs, [train_size, val_size])

    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=shuffle)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)

    context.update({
        "train_loader": train_loader,
        "val_loader": val_loader,
        "y_mean": y_mean,
        "y_std": y_std,
    })

    return {
        "train_loader": train_loader,
        "val_loader": val_loader,
        "y_mean": y_mean,
        "y_std": y_std,
    }

@register_task("train_pose_ranker", category="Pose ranker", description="Train Pose Ranker model with grid search.")
def train_model(config: dict, context: dict):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    train_loader = context.get("train_loader")
    val_loader = context.get("val_loader")
    y_mean = context.get("y_mean")
    y_std = context.get("y_std")
    graphs = context.get("graphs")
    
    if not train_loader or not val_loader or not graphs:
        raise ValueError("Missing train_loader, val_loader, or graphs in context.")

    # Hyperparameters from config with defaults
    param_grid = config.get("param_grid", {
        "hidden_dim": [64, 128, 256],
        "lr": [1e-3, 5e-4, 1e-4]
    })
    num_epochs = config.get("num_epochs", 100)
    early_stop_patience = config.get("early_stop_patience", 5)
    n_outputs = config.get("n_outputs", 9)

    in_dim = graphs[0].x.shape[1]
    results = []

    # training loop begins
    
    for hidden_dim, lr in itertools.product(param_grid["hidden_dim"], param_grid["lr"]):
        logger.info(f"\n=== Grid Search: hidden_dim={hidden_dim}, lr={lr} ===")

        model = PoseScoringGNN(in_dim, hidden_dim=hidden_dim, n_outputs=n_outputs).to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=lr)
        criterion = torch.nn.MSELoss()
        early_stopper = EarlyStopping(patience=early_stop_patience)

        train_losses, val_losses = [], []

        for epoch in range(num_epochs):
            train_loss = train_epoch(model, train_loader, optimizer, criterion)
            val_loss = eval_epoch(model, val_loader, criterion)

            train_losses.append(train_loss)
            val_losses.append(val_loss)

            logger.info(f"Epoch {epoch+1}/{num_epochs} - Train Loss: {train_loss:.4f} - Val Loss: {val_loss:.4f}")

            early_stopper.step(val_loss)
            if early_stopper.should_stop:
                logger.info(f"Early stopping at epoch {epoch+1}")
                break

        results.append({
            "hidden_dim": hidden_dim,
            "lr": lr,
            "best_val_loss": min(val_losses),
            "train_losses": train_losses,
            "val_losses": val_losses,
            "model_state_dict": model.state_dict()
        })

        # Save plot of training curves
        plt.figure(figsize=(6, 4))
        plt.plot(train_losses, label="Train Loss")
        plt.plot(val_losses, label="Validation Loss")
        plt.xlabel("Epoch")
        plt.ylabel("MSE Loss")
        plt.title(f"Loss Curves (hidden_dim={hidden_dim}, lr={lr})")
        plt.legend()
        plt.tight_layout()
        filename = f"_images/loss_curve_hidden{hidden_dim}_lr{lr}.png"
        plt.savefig(filename, dpi=150)
        plt.close()
        logger.debug(f"Saved loss curve to {filename}")

    # Find best result and retrain briefly if desired
    best_result = min(results, key=lambda x: x["best_val_loss"])
    logger.info("\n=== Best Hyperparameters ===")
    logger.info(best_result)

    # Optional final retrain (e.g., 20 epochs)
    best_model = PoseScoringGNN(in_dim, hidden_dim=best_result["hidden_dim"], n_outputs=n_outputs).to(device)
    best_optimizer = torch.optim.Adam(best_model.parameters(), lr=best_result["lr"])
    criterion = torch.nn.MSELoss()

    final_epochs = config.get("final_epochs", 20)
    for epoch in range(final_epochs):
        train_loss = train_epoch(best_model, train_loader, best_optimizer, criterion)
        val_loss = eval_epoch(best_model, val_loader, criterion)
        logger.info(f"Final Training Epoch {epoch+1}/{final_epochs} - Train: {train_loss:.4f}, Val: {val_loss:.4f}")

    model_path = config.get("model_path", "pose_scoring_gnn_best.pth")
    torch.save(best_model.state_dict(), model_path)
    logger.info(f"Saved best model to '{model_path}'")

    # Update context
    context.update({
        "best_model": best_model,
        "train_results": results,
        "best_hyperparams": best_result,
        "model_path": model_path,
    })

    return {
        "best_model": best_model,
        "train_results": results,
        "best_hyperparams": best_result,
        "model_path": model_path,
    }

@register_task("evaluate_model", category="Pose ranker", description="Evaluate trained model stats and plots.")
def evaluate_model(config: dict, context: dict):
    """
    Evaluate the trained model on validation data, print metrics per output, and save parity plots.

    Requires in context:
      - best_model or model
      - val_loader
      - criterion (optional, default MSELoss)
      - device
      - y_mean, y_std for denormalisation

    Saves plots in _images/
    """
    model = context.get("best_model") or context.get("model")
    val_loader = context.get("val_loader")
    criterion = context.get("criterion") or torch.nn.MSELoss()
    device = context.get("device")
    y_mean = context.get("y_mean")
    y_std = context.get("y_std")

    if not model or not val_loader or device is None or y_mean is None or y_std is None:
        raise ValueError("Missing required context entries for evaluation (model, val_loader, device, y_mean, y_std)")

    model.eval()
    total_loss = 0
    all_preds, all_targets = [], []
    with torch.no_grad():
        for batch in val_loader:
            batch = batch.to(device)
            pred = model(batch)
            loss = criterion(pred, batch.y)
            total_loss += loss.item()
            all_preds.append(pred.cpu())
            all_targets.append(batch.y.cpu())

    val_loss = total_loss / len(val_loader)
    all_preds = torch.cat(all_preds).numpy()
    all_targets = torch.cat(all_targets).numpy()

    all_preds = denormalise(all_preds, y_mean, y_std)
    all_targets = denormalise(all_targets, y_mean, y_std)

    metric_names = [
        "ΔRMSD (Å)",
        "Interface RMSD (Å)",
        "f_native",
        "Q-score",
        "ΔSASA",
        "Energy Proxy",
        "H-bonds: count",
        "H-bonds: mean dist",
        "H-bonds: mean angle",
    ]

    n_outputs = all_preds.shape[1]
    if len(metric_names) < n_outputs:
        metric_names += [f"Output_{i+1}" for i in range(len(metric_names), n_outputs)]

    logger.info("\n=== Validation Metrics per Output ===")
    for i, name in enumerate(metric_names[:n_outputs]):
        y_true = all_targets[:, i]
        y_pred = all_preds[:, i]

        mse = mean_squared_error(y_true, y_pred)
        rmse = np.sqrt(mse)
        mae = np.mean(np.abs(y_pred - y_true))
        r2 = r2_score(y_true, y_pred)

        logger.info(f"\n{name}:")
        logger.info(f"  RMSE: {rmse:.4f}")
        logger.info(f"  MAE:  {mae:.4f}")
        logger.info(f"  R²:   {r2:.4f}")

    save_dir = "_images"
    os.makedirs(save_dir, exist_ok=True)

    for i, name in enumerate(metric_names[:n_outputs]):
        y_true = all_targets[:, i]
        y_pred = all_preds[:, i]

        if np.allclose(y_true, y_true[0]):
            logger.info(f"Skipping {name}: constant target values.")
            continue

        reg = LinearRegression().fit(y_true.reshape(-1, 1), y_pred)
        y_fit = reg.predict(y_true.reshape(-1, 1))
        r2_local = r2_score(y_true, y_pred)
        slope, intercept = reg.coef_[0], reg.intercept_

        plt.figure(figsize=(5, 5))
        plt.scatter(y_true, y_pred, alpha=0.6, edgecolor='k', linewidth=0.3)
        plt.plot([y_true.min(), y_true.max()], [y_true.min(), y_true.max()], 'r--', lw=2, label='Ideal')
        plt.plot(y_true, y_fit, 'b-', lw=1.5, label=f'Fit: y = {slope:.2f}x + {intercept:.2f}')
        plt.xlabel(f"True {name}")
        plt.ylabel(f"Predicted {name}")
        plt.title(f"Parity Plot for {name}")
        plt.legend()
        plt.grid(True)
        plt.text(
            0.05, 0.95,
            f"R² = {r2_local:.3f}",
            transform=plt.gca().transAxes,
            fontsize=10,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7)
        )
        plt.tight_layout()

        fname = os.path.join(save_dir, f"parity_{name.replace(' ', '_').replace('(', '').replace(')', '')}.png")
        plt.savefig(fname, dpi=300, bbox_inches='tight')
        plt.close()

    logger.info(f"\nSaved diagnostic plots to: {os.path.abspath(save_dir)}")

    return {"val_loss": val_loss, "predictions": all_preds, "targets": all_targets}
