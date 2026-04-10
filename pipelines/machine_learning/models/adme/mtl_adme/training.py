import itertools
import torch
import numpy as np
import matplotlib.pyplot as plt

from tqdm import trange
from torch_geometric.loader import DataLoader
from torch_geometric.data import Batch
import torch.nn.functional as F
from sklearn.model_selection import train_test_split
from pathlib import Path

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)


# =========================
# COLLATE
# =========================
def custom_collate(batch_list):
    return Batch.from_data_list(batch_list)


# =========================
# LOSSES
# =========================
def masked_mse_loss(pred, target):
    mask = ~torch.isnan(target)

    if mask.sum() == 0:
        return torch.tensor(0.0, device=pred.device, requires_grad=True)

    diff = pred - target
    diff = diff[mask]   # apply mask after subtraction

    return torch.mean(diff ** 2)


def per_task_mse(pred, target):
    mask = ~torch.isnan(target)
    num_tasks = target.shape[1]

    losses = torch.zeros(num_tasks, device=pred.device)
    counts = torch.zeros(num_tasks, device=pred.device)

    for t in range(num_tasks):
        m = mask[:, t]
        if m.sum() > 0:
            losses[t] = torch.sum((pred[m, t] - target[m, t]) ** 2)
            counts[t] = m.sum()

    return losses, counts

def multitask_loss(pred, target, log_vars):
    total_loss = 0.0

    for i in range(pred.shape[1]):
        mask = ~torch.isnan(target[:, i])
        if mask.sum() == 0:
            continue

        pred_i = pred[mask, i]
        target_i = target[mask, i]

        mse = F.mse_loss(pred_i, target_i)

        precision = torch.exp(-log_vars[i])
        total_loss += precision * mse + log_vars[i]

    return total_loss

# =========================
# TRAIN / EVAL EPOCHS
# =========================
def train_epoch(model, loader, optimizer, device):
    model.train()
    total_loss = 0

    num_tasks = None
    task_losses = None
    task_counts = None

    for batch in loader:
        batch = batch.to(device)

        if hasattr(batch, "global_features"):
            gf = batch.global_features.to(device).float()
            if gf.dim() == 1:
                gf = gf.unsqueeze(0).repeat(batch.num_graphs, 1)
            batch.global_features = gf

        optimizer.zero_grad()

        out = model(batch)
        target = batch.y.float()

        loss = multitask_loss(out, target, model.log_vars)
        loss.backward()
        optimizer.step()

        total_loss += loss.item()

        # per-task tracking
        if num_tasks is None:
            num_tasks = target.shape[1]
            task_losses = torch.zeros(num_tasks, device=device)
            task_counts = torch.zeros(num_tasks, device=device)

        losses, counts = per_task_mse(out, target)
        task_losses += losses
        task_counts += counts

    task_mse = (task_losses / (task_counts + 1e-8)).detach().cpu().numpy()

    return total_loss / len(loader), task_mse


def eval_epoch(model, loader, device):
    model.eval()
    total_loss = 0

    all_preds, all_labels = [], []

    num_tasks = None
    task_losses = None
    task_counts = None

    with torch.no_grad():
        for batch in loader:
            batch = batch.to(device)

            if hasattr(batch, "global_features"):
                gf = batch.global_features.to(device).float()
                if gf.dim() == 1:
                    gf = gf.unsqueeze(0).repeat(batch.num_graphs, 1)
                batch.global_features = gf

            out = model(batch)
            target = batch.y.float()

            loss = multitask_loss(out, target, model.log_vars)
            total_loss += loss.item()

            all_preds.append(out.cpu())
            all_labels.append(batch.y.cpu())

            if num_tasks is None:
                num_tasks = target.shape[1]
                task_losses = torch.zeros(num_tasks, device=device)
                task_counts = torch.zeros(num_tasks, device=device)

            losses, counts = per_task_mse(out, target)
            task_losses += losses
            task_counts += counts

    task_mse = (task_losses / (task_counts + 1e-8)).detach().cpu().numpy()

    return (
        total_loss / len(loader),
        torch.cat(all_preds),
        torch.cat(all_labels),
        task_mse
    )


# =========================
# MAIN TRAIN FUNCTION (PIPELINE ENTRY)
# =========================
def train(context, config):
    data_list = context["graphs"]
    device = context["device"]
    ModelClass = context["model_class"]

    logger.debug(f"Number of graphs: {len(data_list)}")

    # -------------------------
    # SPLIT
    # -------------------------
    train_list, val_list = train_test_split(
        data_list,
        test_size=config.get("val_fraction", 0.2),
        random_state=config.get("random_seed", 42),
    )

    train_loader = DataLoader(
        train_list,
        batch_size=config["batch_size"],
        shuffle=True,
        collate_fn=custom_collate,
    )

    val_loader = DataLoader(
        val_list,
        batch_size=config["batch_size"],
        shuffle=False,
        collate_fn=custom_collate,
    )

    # -------------------------
    # DIMENSIONS
    # -------------------------
    sample = data_list[0]

    input_dim = sample.x.shape[1]
    edge_dim = sample.edge_attr.shape[1]
    global_feat_dim = sample.global_features.shape[-1]
    fp_dim = sample.fp.shape[-1]
    num_tasks = sample.y.shape[-1]

    # -------------------------
    # OUTPUT DIRS
    # -------------------------
    loss_curve_dir = Path(config.get("loss_curve_dir", "outputs/mtl/loss_curves"))
    loss_curve_dir.mkdir(parents=True, exist_ok=True)

    model_path = Path(config.get("model_path", "outputs/mtl/model.pth"))
    model_path.parent.mkdir(parents=True, exist_ok=True)

    best_model, best_loss, best_params = None, float("inf"), None

    # =========================
    # GRID SEARCH
    # =========================
    for lr, hidden_dim in itertools.product(
        config["param_grid"]["lr"],
        config["param_grid"]["hidden_dim"],
    ):
        logger.info(f"\nTraining with lr={lr}, hidden_dim={hidden_dim}")

        model = ModelClass(
            input_dim=input_dim,
            edge_dim=edge_dim,
            global_feat_dim=global_feat_dim,
            fp_dim=fp_dim,
            hidden_dim=hidden_dim,
            num_tasks=num_tasks,
        ).to(device)

        optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-5)

        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, mode="min", factor=0.5, patience=10, verbose=True
        )

        best_val_loss = float("inf")
        patience = config.get("patience", 50)
        counter = 0

        train_losses, val_losses = [], []
        train_task_hist, val_task_hist = [], []

        best_state = {k: v.detach().cpu() for k, v in model.state_dict().items()}

        # =========================
        # TRAINING LOOP
        # =========================
        for epoch in trange(
            config.get("max_epochs", 200),
            desc=f"lr={lr}, hidden={hidden_dim}",
        ):
            train_loss, train_task = train_epoch(model, train_loader, optimizer, device)

            val_loss, _, _, val_task = eval_epoch(model, val_loader, device)

            train_losses.append(train_loss)
            val_losses.append(val_loss)

            train_task_hist.append(train_task)
            val_task_hist.append(val_task)

            scheduler.step(val_loss)

            if not np.isnan(val_loss) and val_loss < best_val_loss:
                best_val_loss = val_loss
                best_state = {k: v.detach().cpu() for k, v in model.state_dict().items()}
                counter = 0
            else:
                counter += 1

            if counter >= patience:
                logger.info("Early stopping triggered")
                break

        # -------------------------
        # LOAD BEST STATE
        # -------------------------
        model.load_state_dict(best_state)
        model.to(device)

        # -------------------------
        # PLOTTING
        # -------------------------
        train_task_hist = np.array(train_task_hist)
        val_task_hist = np.array(val_task_hist)

        # Overall loss
        plt.figure()
        plt.plot(train_losses, label="Train")
        plt.plot(val_losses, label="Val")
        plt.legend()
        plt.title(f"lr={lr}, hidden={hidden_dim}")
        plt.savefig(loss_curve_dir / f"overall_lr{lr}_hd{hidden_dim}.png", dpi=150)
        plt.close()

        # Per-task loss
        n_tasks = num_tasks
        n_cols = 3
        n_rows = int(np.ceil(n_tasks / n_cols))

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
        axes = axes.flatten()

        for i in range(n_tasks):
            axes[i].plot(train_task_hist[:, i], label="Train")
            axes[i].plot(val_task_hist[:, i], label="Val")
            axes[i].set_title(f"Task {i}")
            axes[i].legend()

        for j in range(i + 1, len(axes)):
            fig.delaxes(axes[j])

        plt.tight_layout()
        plt.savefig(loss_curve_dir / f"per_task_lr{lr}_hd{hidden_dim}.png", dpi=150)
        plt.close()

        # -------------------------
        # FINAL EVAL
        # -------------------------
        final_val_loss, val_preds, val_labels, _ = eval_epoch(model, val_loader, device)

        y_pred = val_preds.cpu().numpy()
        y_true = val_labels.cpu().numpy()

        pred_dir = Path(config.get("predictions_dir", "outputs/mtl_adme/preds"))
        pred_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Saved predictions to {pred_dir}")

        if np.isnan(final_val_loss):
            logger.warning("Skipping model due to NaN val_loss")
            continue

        if final_val_loss < best_loss:
            best_loss = final_val_loss
            best_model = model
            best_params = {"lr": lr, "hidden_dim": hidden_dim}

            torch.save(
                {
                    "model_state_dict": model.state_dict(),
                    "input_dim": input_dim,
                    "edge_dim": edge_dim,
                    "global_feat_dim": global_feat_dim,
                    "fp_dim": fp_dim,
                    "num_tasks": num_tasks,
                },
                model_path,
            )

            logger.info(f"Saved best model (val_loss={best_loss:.4f})")

            pred_dir = Path(config.get("predictions_dir", "outputs/mtl_adme/preds"))
            pred_dir.mkdir(parents=True, exist_ok=True)

            np.save(pred_dir / "y_pred.npy", y_pred)
            np.save(pred_dir / "y_true.npy", y_true)

            logger.info(f"Saved best predictions to {pred_dir}")

    logger.info(f"\nBest params: {best_params}")

    if best_model is None:
        raise RuntimeError("Training failed: no valid model found (all val losses NaN)")

    return {
        "y_pred": y_pred,
        "y_true": y_true,
        "task_names": context.get("task_names"),
        "predictions_dir": str(pred_dir),
        "best_val_loss": best_loss,
        "best_params": best_params,
        "model_path": str(model_path),
    }