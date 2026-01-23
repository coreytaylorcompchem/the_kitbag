import itertools
import torch
import matplotlib.pyplot as plt
from tqdm import trange
from torch_geometric.loader import DataLoader
from torch_geometric.data import Batch
from sklearn.model_selection import train_test_split
from pathlib import Path

from models.adme.cyp3a4.model import GINRegressor

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)


# Collate must preserve MoleculeData
def custom_collate(batch_list):
    return Batch.from_data_list(batch_list)


def train_epoch(model, loader, optimizer, loss_fn, device):
    model.train()
    total_loss = 0.0
    for batch in loader:
        batch = batch.to(device)

        # Move global_features to device and ensure correct shape
        if hasattr(batch, 'global_features'):
            global_feat = batch.global_features.to(device).float()
            if global_feat.dim() == 1:
                global_feat = global_feat.unsqueeze(0).repeat(batch.num_graphs, 1)
            batch.global_features = global_feat

        # Ensure labels are float32 and NaN-free
        batch.y = torch.nan_to_num(batch.y, nan=0.0).float()

        optimizer.zero_grad()
        out = model(batch)
        loss = loss_fn(out, batch.y)
        loss.backward()
        optimizer.step()
        total_loss += loss.item()

    return total_loss / len(loader)


def eval_epoch(model, loader, loss_fn, device):
    model.eval()
    total_loss = 0.0
    all_preds, all_labels = [], []

    with torch.no_grad():
        for batch in loader:
            batch = batch.to(device)

            if hasattr(batch, 'global_features'):
                global_feat = batch.global_features.to(device).float()
                if global_feat.dim() == 1:
                    global_feat = global_feat.unsqueeze(0).repeat(batch.num_graphs, 1)
                batch.global_features = global_feat

            batch.y = torch.nan_to_num(batch.y, nan=0.0).float()

            out = model(batch)
            loss = loss_fn(out, batch.y)
            total_loss += loss.item()
            all_preds.append(out.cpu())
            all_labels.append(batch.y.view(-1).cpu())

    return total_loss / len(loader), torch.cat(all_preds), torch.cat(all_labels)


def train(context, config):
    """
    Pipeline-compatible training using MoleculeData.
    Uses global_feat_dim already set in context by featurisation.
    Includes hyperparameter search, early stopping, loss curve plotting, and model saving.
    """

    data_list = context["graphs"]
    device = context["device"]

    global_feat_dim = context.get("global_feat_dim", 0)

    # Split data
    train_list, val_list = train_test_split(
        data_list,
        test_size=config.get("val_fraction", 0.2),
        random_state=config.get("random_seed", 42),
    )

    train_loader = DataLoader(
        train_list,
        batch_size=config["batch_size"],
        shuffle=True,
        collate_fn=custom_collate
    )
    val_loader = DataLoader(
        val_list,
        batch_size=config["batch_size"],
        shuffle=False,
        collate_fn=custom_collate
    )

    input_dim = data_list[0].x.shape[1]

    best_model, best_loss, best_params = None, float("inf"), None
    loss_curve_dir = Path(config.get("loss_curve_dir", "_plots"))
    loss_curve_dir.mkdir(parents=True, exist_ok=True)

    # Prepare model output path
    model_path = config.get("model_path", "outputs/models/best_model.pth")
    model_path = Path(model_path)
    model_path.parent.mkdir(parents=True, exist_ok=True)

    for lr, hidden_dim in itertools.product(
        config["param_grid"]["lr"], config["param_grid"]["hidden_dim"]
    ):
        model = GINRegressor(
            input_dim=input_dim,
            hidden_dim=hidden_dim,
            global_feat_dim=global_feat_dim
        ).to(device)

        optimizer = torch.optim.Adam(model.parameters(), lr=lr)
        loss_fn = torch.nn.MSELoss()

        best_val_loss = float("inf")
        patience, counter = config.get("patience", 10), 0
        train_losses, val_losses = [], []

        # Save initial best state
        best_state = {k: v.detach().cpu() for k, v in model.state_dict().items()}

        for epoch in trange(config.get("max_epochs", 100), desc=f"Training lr={lr}, hidden={hidden_dim}"):
            train_loss = train_epoch(model, train_loader, optimizer, loss_fn, device)
            val_loss, _, _ = eval_epoch(model, val_loader, loss_fn, device)

            train_losses.append(train_loss)
            val_losses.append(val_loss)
            logger.info(f"Epoch {epoch+1} | Train: {train_loss:.4f} | Val: {val_loss:.4f}")

            # Early stopping
            if val_loss < best_val_loss:
                best_val_loss = val_loss
                best_state = {k: v.detach().cpu() for k, v in model.state_dict().items()}
                counter = 0
            else:
                counter += 1

            if counter >= patience:
                logger.info("Early stopping triggered.")
                break

        # Restore best state
        model.load_state_dict(best_state)
        model.to(device)

        # Save loss curves
        plt.figure(figsize=(8, 5))
        plt.plot(range(1, len(train_losses) + 1), train_losses, label="Train Loss")
        plt.plot(range(1, len(val_losses) + 1), val_losses, label="Validation Loss")
        plt.xlabel("Epoch")
        plt.ylabel("MSE Loss")
        plt.title(f"Loss Curves lr={lr}, hidden={hidden_dim}")
        plt.legend()
        plt.tight_layout()
        plt.savefig(loss_curve_dir / f"loss_curve_lr{lr}_hidden{hidden_dim}.png", dpi=150)
        plt.close()

        # Update global best
        final_val_loss, _, _ = eval_epoch(model, val_loader, loss_fn, device)
        if final_val_loss < best_loss:
            best_loss = final_val_loss
            best_model = model
            best_params = {"lr": lr, "hidden_dim": hidden_dim}

            # Save the best model to disk
            torch.save({
                "model_state_dict": model.state_dict(),
                "input_dim": input_dim,
                "hidden_dim": hidden_dim,
                "global_feat_dim": global_feat_dim
                }, model_path)
            
            logger.info(f"Saved best model (val_loss={best_loss:.4f}) to {model_path}")

    logger.info(f"\nBest params: {best_params}")
    return {
        "model": best_model,
        "best_val_loss": best_loss,
        "best_params": best_params,
        "model_path": str(model_path)
    }

