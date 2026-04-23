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


########## HELPERS ##########

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

def multitask_loss(pred, target, log_vars, active_tasks=None):
    total_loss = None

    for i in range(pred.shape[1]):

        if active_tasks is not None and i not in active_tasks:
            continue

        mask = ~torch.isnan(target[:, i])
        if mask.sum() == 0:
            continue

        pred_i = pred[mask, i]
        target_i = target[mask, i]

        mse = F.mse_loss(pred_i, target_i)

        log_var = torch.clamp(log_vars[i], -5.0, 5.0)
        precision = torch.exp(-log_var)

        loss_i = precision * mse + log_var

        if total_loss is None:
            total_loss = loss_i
        else:
            total_loss = total_loss + loss_i

    # fallback
    if total_loss is None:
        return pred.sum() * 0.0

    return total_loss

# def compute_task_losses(pred, target, log_vars, active_tasks=None):
#     task_losses = []

#     for i in range(pred.shape[1]):

#         if active_tasks is not None and i not in active_tasks:
#             continue

#         mask = ~torch.isnan(target[:, i])
#         if mask.sum() == 0:
#             continue

#         pred_i = pred[mask, i]
#         target_i = target[mask, i]

#         mse = F.mse_loss(pred_i, target_i)

#         log_var = torch.clamp(log_vars[i], -5.0, 5.0)
#         precision = torch.exp(-log_var)

#         loss_i = precision * mse + log_var
#         task_losses.append(loss_i)

#     return task_losses

########### PCGrad implmentation

class PCGrad:
    def __init__(self, optimizer):
        self.optimizer = optimizer

    def pc_backward(self, model, out, target, loss_fns):
        """
        Memory-efficient PCGrad:
        - single forward pass
        - multiple backward passes with retain_graph
        """

        params = [p for p in model.parameters() if p.requires_grad]
        final_grads = [torch.zeros_like(p) for p in params]

        for i, loss_fn in enumerate(loss_fns):

            self.optimizer.zero_grad(set_to_none=True)

            loss_i = loss_fn(out, target)
            if loss_i is None:
                continue

            retain = i < (len(loss_fns) - 1)

            loss_i.backward(retain_graph=retain)

            grad_i = [
                p.grad.clone() if p.grad is not None else None
                for p in params
            ]

            for k in range(len(final_grads)):
                if grad_i[k] is not None:
                    final_grads[k] += grad_i[k]

        self.optimizer.zero_grad(set_to_none=True)

        for p, g in zip(params, final_grads):
            if g is not None:
                p.grad = g / len(loss_fns)

    def step(self):
        self.optimizer.step()

# =========================
# TRAIN / EVAL EPOCHS
# =========================
def train_epoch(model, loader, optimizer, device, active_tasks=None):

    all_preds, all_labels = [], []

    model.train()
    total_loss = 0

    num_tasks = None
    task_losses = None
    task_counts = None

    pcgrad = PCGrad(optimizer)

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

        use_pcgrad = getattr(model, "use_pcgrad", True)

        if use_pcgrad:

            # -------------------------
            # Select task indices
            # -------------------------
            task_indices = [
                i for i in range(out.shape[1])
                if (active_tasks is None or i in active_tasks)
            ]

            max_tasks = getattr(model, "pcgrad_max_tasks", 2)
            task_indices = task_indices[:max_tasks]

            # -------------------------
            # Build loss functions
            # -------------------------
            def make_loss_fn(task_idx):
                def loss_fn(out, target):
                    mask = ~torch.isnan(target[:, task_idx])
                    if mask.sum() == 0:
                        return None

                    pred_i = out[mask, task_idx]
                    target_i = target[mask, task_idx]

                    mse = F.mse_loss(pred_i, target_i)

                    log_var = torch.clamp(model.log_vars[task_idx], -5.0, 5.0)
                    precision = torch.exp(-log_var)

                    return precision * mse + log_var
                return loss_fn

            loss_fns = [make_loss_fn(i) for i in task_indices]

            if len(loss_fns) > 0:
                # --- PCGrad update ---
                pcgrad.pc_backward(model, out, target, loss_fns)
                pcgrad.step()

                # --- REAL LOSS FOR LOGGING ---
                with torch.no_grad():
                    losses = []
                    for loss_fn in loss_fns:
                        l = loss_fn(out, target)
                        if l is not None:
                            losses.append(l)

                    if len(losses) > 0:
                        loss = torch.stack(losses).mean()
                    else:
                        loss = torch.tensor(0.0, device=device)
            else:
                loss = torch.tensor(0.0, device=device)
        else:
            loss = multitask_loss(out, target, model.log_vars, active_tasks)

            if loss.requires_grad:
                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                optimizer.step()
            else:
                # No valid tasks in this batch; skip update safely
                optimizer.zero_grad(set_to_none=True)

        total_loss += loss.item()

        # per-task tracking
        if num_tasks is None:
            num_tasks = target.shape[1]
            task_losses = torch.zeros(num_tasks, device=device)
            task_counts = torch.zeros(num_tasks, device=device)

        losses, counts = per_task_mse(out, target)
        task_losses += losses
        task_counts += counts

        all_preds.append(out.detach())
        all_labels.append(target.detach())

    task_mse = (task_losses / (task_counts + 1e-8))
    task_rmse = torch.sqrt(task_mse).detach().cpu().numpy()

    task_r2 = per_task_r2(torch.cat(all_preds), torch.cat(all_labels)).detach().cpu().numpy()

    return (
        total_loss / len(loader),
        torch.cat(all_preds),
        torch.cat(all_labels),
        task_rmse,
        task_r2
    )


def eval_epoch(model, loader, device, active_tasks=None):
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

            loss = multitask_loss(out, target, model.log_vars, active_tasks)
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

def set_requires_grad(module, requires_grad):
    for p in module.parameters():
        p.requires_grad = requires_grad

def apply_curriculum_freezing(model, stage_idx, curriculum_cfg, context):
    strategy = curriculum_cfg.get("strategy", "none")

    if strategy == "none":
        return

    # First freeze everything
    for p in model.parameters():
        p.requires_grad = False

    # Always allow shared trunk + backbone
    set_requires_grad(model.input_proj, True)
    set_requires_grad(model.convs, True)
    set_requires_grad(model.norms, True)
    set_requires_grad(model.shared, True)
    set_requires_grad(model.fp_mlp, True)
    set_requires_grad(model.global_proj, True)

    stages = curriculum_cfg["stages"]

    if strategy == "freeze_previous":
        active_stages = [stage_idx]

    elif strategy == "progressive_unfreeze":
        active_stages = list(range(stage_idx + 1))

    else:
        raise ValueError(f"Unknown curriculum strategy: {strategy}")

    # Enable group trunks + heads for active stages
    for i in active_stages:
        stage = stages[i]
        group_name = stage["name"]
        if isinstance(stage["tasks"][0], str):
            task_indices = [
                context["task_names"].index(t)
                for t in stage["tasks"]
            ]
        else:
            task_indices = stage["tasks"]

        set_requires_grad(model.group_trunks[group_name], True)

        for t in task_indices:
            set_requires_grad(model.heads[t], True)

def filter_dataset_for_tasks(data_list, task_indices):
    filtered = []

    for data in data_list:
        y = data.y

        # Ensure shape is [num_tasks]
        if y.dim() == 2:
            y = y.squeeze(0)

        if not torch.isnan(y[task_indices]).all():
            filtered.append(data)

    return filtered

def per_task_r2(pred, target):
    mask = ~torch.isnan(target)
    num_tasks = target.shape[1]

    r2 = torch.full((num_tasks,), float("nan"), device=pred.device)

    for t in range(num_tasks):
        m = mask[:, t]
        if m.sum() > 1:
            y_true = target[m, t]
            y_pred = pred[m, t]

            ss_res = torch.sum((y_true - y_pred) ** 2)
            ss_tot = torch.sum((y_true - torch.mean(y_true)) ** 2)

            if ss_tot > 0:
                r2[t] = 1 - ss_res / ss_tot

    return r2

######## TRAINING LOOPS ##########

def train_curriculum(context, config, params):

    device = context["device"]
    ModelClass = context["model_class"]
    data_list = context["graphs"]

    curriculum_cfg = config["curriculum"]
    stages = curriculum_cfg["stages"]
    strategy = curriculum_cfg.get("strategy", "none")

    # Dataloaders
    train_list, val_list = train_test_split(
        data_list,
        test_size=config.get("val_fraction", 0.2),
        random_state=config.get("random_seed", 42),
    )

   # Model init
    sample = data_list[0]

    model = ModelClass(
        input_dim=sample.x.shape[1],
        edge_dim=sample.edge_attr.shape[1],
        global_feat_dim=sample.global_features.shape[-1],
        fp_dim=sample.fp.shape[-1],
        num_tasks=sample.y.shape[-1],
        task_groups=context.get("task_groups"),
        **{k: v for k, v in params.items() if k != "lr"}
    ).to(device)

    model.use_pcgrad = config.get("use_pcgrad", True)

    if model.use_pcgrad:
        logger.info("PCGrad ENABLED for training")
    else:
        logger.info("PCGrad DISABLED (using standard multitask loss)")

    best_loss = float("inf")
    best_model_state = None 

    logger.info(f"Curriculum freezing strategy applied: {strategy}")

    for stage_idx, stage in enumerate(stages):

        logger.debug(f"[Stage: {stage['name']}] Available group trunks: {list(model.group_trunks.keys())}")

        if isinstance(stage["tasks"][0], str):
            stage_tasks = [
                context["task_names"].index(t)
                for t in stage["tasks"]
            ]
        else:
            stage_tasks = stage["tasks"]
        stage_epochs = stage["epochs"]

         # Stage-specific dataset
        stage_train_list = filter_dataset_for_tasks(train_list, stage_tasks)
        stage_val_list = filter_dataset_for_tasks(val_list, stage_tasks)
        # stage_val_list = val_list

        stage_train_list = train_list
        stage_val_list = val_list

        if len(stage_train_list) == 0:
            logger.warning(f"No data for stage {stage['name']}, skipping")
            continue

        train_loader = DataLoader(
            stage_train_list,
            batch_size=config["batch_size"],
            shuffle=True,
            collate_fn=custom_collate,
        )

        val_loader = DataLoader(
            stage_val_list,
            batch_size=config["batch_size"],
            shuffle=False,
            collate_fn=custom_collate,
        )

        # Stage-specific params
        lr = params.get("lr", 1e-3)

        optimizer = torch.optim.Adam(
            filter(lambda p: p.requires_grad, model.parameters()),
            lr=lr
            )
        patience = stage.get("patience", 10)

        logger.info(f"=== Curriculum stage: {stage['name']} | tasks={stage_tasks} ===")

        # History tracking
        train_losses = []
        val_losses = []
        train_task_hist = []
        val_task_hist = []

        apply_curriculum_freezing(model, stage_idx, curriculum_cfg, context)

        if stage_idx == 0 or curriculum_cfg.get("reset_optimizer", True):
            optimizer = torch.optim.Adam(
                filter(lambda p: p.requires_grad, model.parameters()),
                lr=lr
            )

        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, mode="min", factor=0.5, patience=5, verbose=True
        )

        # Per-stage early stopping
        best_stage_loss = float("inf")
        best_stage_state = None
        counter = 0

        for epoch in range(stage_epochs):

            # =========================
            # TRAIN
            # =========================
            train_loss, train_preds, train_labels, train_rmse, train_r2 = train_epoch(
                model, train_loader, optimizer, device, active_tasks=stage_tasks
            )

            # =========================
            # MASK TRAIN RMSE TO STAGE TASKS ONLY
            # =========================
            train_rmse_stage = train_rmse.copy()
            mask_arr = np.ones_like(train_rmse_stage, dtype=bool)
            mask_arr[stage_tasks] = False
            train_rmse_stage[mask_arr] = np.nan

            # =========================
            # VALIDATION
            # =========================
            val_loss, val_preds, val_labels, val_mse = eval_epoch(
                model, val_loader, device, active_tasks=stage_tasks
            )

            val_rmse = np.sqrt(val_mse)

            val_rmse_stage = val_rmse.copy()
            mask_arr = np.ones_like(val_rmse_stage, dtype=bool)
            mask_arr[stage_tasks] = False
            val_rmse_stage[mask_arr] = np.nan

            val_r2 = per_task_r2(
                val_preds.to(device),
                val_labels.to(device)
            ).detach().cpu().numpy()

            # =========================
            # STORE HISTORY (NEW)
            # =========================
            train_losses.append(train_loss)
            val_losses.append(val_loss)
            train_task_hist.append(train_rmse_stage)
            val_task_hist.append(val_rmse_stage)

            logger.info(
                f"Stage {stage['name']} | Epoch {epoch} | "
                f"Train Loss {train_loss:.4f} | Val Loss {val_loss:.4f}"
            )

            logger.debug(f"Val RMSE: {np.round(val_rmse_stage, 3)}")
            logger.debug(f"Val R2:   {np.round(val_r2, 3)}")

            # =========================
            # EARLY STOPPING
            # =========================
            if val_loss < best_stage_loss:
                best_stage_loss = val_loss
                best_stage_state = {
                    k: v.detach().cpu() for k, v in model.state_dict().items()
                }
                counter = 0
            else:
                counter += 1

            if counter >= patience:
                logger.info(f"Early stopping in stage {stage['name']}")
                break

        # Restore best weights for this stage
        if best_stage_state is not None:
            model.load_state_dict(best_stage_state)
            model.to(device)
        
        # =========================
        # PLOTTING
        # =========================
        train_task_hist = np.array(train_task_hist)
        val_task_hist = np.array(val_task_hist)

        loss_curve_dir = Path(config.get("loss_curve_dir", "outputs/mtl_curriculum/loss_curves"))
        loss_curve_dir.mkdir(parents=True, exist_ok=True)

        stage_name = stage["name"]

        # Overall loss
        plt.figure()
        plt.plot(train_losses, label="Train")
        plt.plot(val_losses, label="Val")
        plt.legend()
        plt.title(f"Stage: {stage_name}")
        param_str = "_".join(f"{k}={v}" for k,v in params.items())
        plt.savefig(loss_curve_dir / f"{param_str}_overall_stage_{stage_name}.png")
        plt.close()

        # Per-task
        n_tasks = train_task_hist.shape[1]
        n_cols = 3
        n_rows = int(np.ceil(n_tasks / n_cols))

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
        axes = axes.flatten()

        task_names = context.get("task_names", [f"task_{i}" for i in range(n_tasks)])

        for i in range(n_tasks):
            axes[i].plot(train_task_hist[:, i], label="Train")
            axes[i].plot(val_task_hist[:, i], label="Val")
            axes[i].set_title(task_names[i])
            axes[i].legend()

        for j in range(i + 1, len(axes)):
            fig.delaxes(axes[j])

        plt.tight_layout()
        plt.savefig(loss_curve_dir / f"per_task_stage_{stage_name}.png", dpi=150)
        plt.close()

        # validation per stage (for saving preds + global best tracking)
        val_loss, val_preds, val_labels, _ = eval_epoch(model, val_loader, device)

        y_pred = val_preds.cpu().numpy()
        y_true = val_labels.cpu().numpy()

        # Save to context
        context["y_pred"] = y_pred
        context["y_true"] = y_true

        # Save to disk
        pred_dir = Path(config.get("predictions_dir", "outputs/mtl_adme/preds"))
        pred_dir.mkdir(parents=True, exist_ok=True)

        np.save(pred_dir / "y_pred.npy", y_pred)
        np.save(pred_dir / "y_true.npy", y_true)

        logger.info(f"Saved predictions to {pred_dir}")

        model_path = Path(config.get("model_path", "outputs/mtl_adme/models/model.pth"))
        model_path.parent.mkdir(parents=True, exist_ok=True)

        torch.save({
            "model_state_dict": best_model_state,
            "params": params,
            "task_names": context.get("task_names"),
        }, model_path)

        logger.info(f"Saved best model to {model_path}")

        # Store best model
        if val_loss < best_loss:
            best_loss = val_loss
            best_model_state = {
                k: v.detach().cpu() for k, v in model.state_dict().items()
            }

    # Restore global best model
    if best_model_state is not None:
        model.load_state_dict(best_model_state)
        model.to(device)

    return {
        "model": model,
        "params": params,
        "best_val_loss": best_loss,
        "y_pred": y_pred,
        "y_true": y_true,
        "task_names": context.get("task_names"),
        "predictions_dir": str(pred_dir),
    }

# =========================
# MAIN TRAIN FUNCTION (PIPELINE ENTRY)
# =========================

def train(context, config):
    training_mode = config.get("training_mode", "multitask")

    if training_mode == "curriculum":

        data_list = context["graphs"]
        sample = data_list[0]
        num_tasks = sample.y.shape[-1]

        task_names = context.get("task_names", [f"task_{i}" for i in range(num_tasks)])

        task_groups_cfg = config.get("task_groups", {})

        if not task_groups_cfg:
            raise ValueError("task_groups must be provided in config for curriculum mode")

        task_groups = {}
        for group, names in task_groups_cfg.items():
            indices = [task_names.index(n) for n in names]
            task_groups[group] = indices

        context["task_groups"] = task_groups
        context["task_names"] = task_names

        logger.info(f"[Curriculum] Resolved task groups: {task_groups}")

        param_grid = config["param_grid"]
        keys = list(param_grid.keys())

        best_result = None
        best_loss = float("inf")

        for values in itertools.product(*param_grid.values()):
            params = dict(zip(keys, values))

            logger.info(f"\n[Curriculum] Training with params: {params}")

            result = train_curriculum(context, config, params)

            if result["best_val_loss"] < best_loss:
                best_loss = result["best_val_loss"]
                best_result = result
                best_result["best_params"] = params

        return best_result
    data_list = context["graphs"]
    device = context["device"]
    ModelClass = context["model_class"]

    logger.debug(f"Number of graphs: {len(data_list)}")

    def format_params(param_dict):
        return "_".join(f"{k}={v}" for k, v in param_dict.items())

    # -------------------------
    # SPLIT # using random splits
    # -------------------------
    train_list, val_list = train_test_split(
        data_list,
        test_size=config.get("val_fraction", 0.2),
        random_state=config.get("random_seed", 42),
    )

    # Use scaffold splits TODO: add this choice to yaml

    # train_list = context["train_loader"].dataset
    # val_list = context["val_loader"].dataset

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

    task_names = context.get("task_names", [f"task_{i}" for i in range(num_tasks)])

    task_groups_cfg = config.get("task_groups", {})

    task_groups = {}
    for group, names in task_groups_cfg.items():
        indices = [task_names.index(n) for n in names]
        task_groups[group] = indices

    context["task_groups"] = task_groups

    logger.info(f"Resolved task groups: {task_groups}")

    # OUTPUT DIRS
    loss_curve_dir = Path(config.get("loss_curve_dir", "outputs/mtl/loss_curves"))
    loss_curve_dir.mkdir(parents=True, exist_ok=True)

    model_path = Path(config.get("model_path", "outputs/mtl/model.pth"))
    model_path.parent.mkdir(parents=True, exist_ok=True)

    best_model, best_loss, best_params = None, float("inf"), None

    # =========================
    # GRID SEARCH
    # =========================

    param_grid = config["param_grid"]
    keys = list(param_grid.keys())

    for values in itertools.product(*param_grid.values()):
        params = dict(zip(keys, values))

        logger.info(f"\nTraining with params: {params}")

        lr = params.get("lr")

        model = ModelClass(
            input_dim=input_dim,
            edge_dim=edge_dim,
            global_feat_dim=global_feat_dim,
            fp_dim=fp_dim,
            num_tasks=num_tasks,
            task_groups=context.get("task_groups"),
            **{k: v for k, v in params.items() if k != "lr"}).to(device)
        
        model.use_pcgrad = config.get("use_pcgrad", True)

        if model.use_pcgrad:
            logger.info("PCGrad ENABLED for training")
        else:
            logger.info("PCGrad DISABLED (using standard multitask loss)")

        optimizer = torch.optim.Adam(
            model.parameters(),
            lr=lr,
            weight_decay=params.get("weight_decay", 1e-5),
        )

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
            desc=str(params),
        ):
            train_loss, _, _, train_task_rmse, train_task_r2 = train_epoch(model, train_loader, optimizer, device)

            val_loss, val_preds, val_labels, val_mse = eval_epoch(model, val_loader, device)
            val_task_rmse = np.sqrt(val_mse)
            scheduler.step(val_loss)

            train_losses.append(train_loss)
            val_losses.append(val_loss)

            train_task_hist.append(train_task_rmse)
            val_task_hist.append(val_task_rmse)

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

        param_str = format_params(params)

        # Overall loss
        plt.figure()
        plt.plot(train_losses, label="Train")
        plt.plot(val_losses, label="Val")
        plt.legend()
        plt.title(param_str)
        plt.savefig(loss_curve_dir / f"overall_{param_str}.png", dpi=150)
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
            axes[i].set_title(task_names[i])
            axes[i].legend()

        for j in range(i + 1, len(axes)):
            fig.delaxes(axes[j])

        plt.tight_layout()
        plt.savefig(loss_curve_dir / f"per_task_{param_str}.png", dpi=150)
        plt.close()

        # -------------------------
        # FINAL EVAL
        # -------------------------
        final_val_loss, val_preds, val_labels, _ = eval_epoch(model, val_loader, device)

        y_pred = val_preds.cpu().numpy()
        y_true = val_labels.cpu().numpy()

        context["y_pred"] = y_pred
        context["y_true"] = y_true

        pred_dir = Path(config.get("predictions_dir", "outputs/mtl_adme/preds"))
        pred_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Saved predictions to {pred_dir}")

        if np.isnan(final_val_loss):
            logger.warning("Skipping model due to NaN val_loss")
            continue

        if final_val_loss < best_loss:
            best_loss = final_val_loss
            best_model = model
            best_params = params.copy()

            torch.save(
                {
                    "model_state_dict": model.state_dict(),
                    "input_dim": input_dim,
                    "edge_dim": edge_dim,
                    "global_feat_dim": global_feat_dim,
                    "fp_dim": fp_dim,
                    "num_tasks": num_tasks,
                    "params": params,
                    "task_names": task_names,
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