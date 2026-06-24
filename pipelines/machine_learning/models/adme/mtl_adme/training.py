import itertools
import torch
import joblib
import json

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

def compute_task_weights(data_list, num_tasks, power=0.5, eps=1e-8):
    """
    Inverse-frequency task weighting.

    More sparse tasks receive larger weights.
    """

    counts = torch.zeros(num_tasks)

    for data in data_list:
        y = data.y

        if y.dim() == 2:
            y = y.squeeze(0)

        counts += (~torch.isnan(y)).float()

    weights = 1.0 / torch.pow(counts + eps, power)

    # normalise mean weight = 1
    weights = weights / weights.mean()

    return weights

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

def multitask_loss(pred, target, active_tasks=None, task_weights=None):

    total_loss = torch.tensor(0.0, device=pred.device)
    count = 0

    for i in range(pred.shape[1]):

        # =====================================
        # CURRICULUM FILTERING
        # =====================================
        if active_tasks is not None and i not in active_tasks:
            continue

        mask = ~torch.isnan(target[:, i])

        if mask.sum() == 0:
            continue

        pred_i = pred[mask, i]
        target_i = target[mask, i]

        loss_i = F.huber_loss(
            pred_i,
            target_i,
            delta=1.0
        )

        if task_weights is not None:
            loss_i = loss_i * task_weights[i]

        total_loss += loss_i
        count += 1

    if count == 0:
        return torch.tensor(0.0, device=pred.device)

    return total_loss / count

# =========================
# LOSS FACTORY
# =========================
def compute_base_loss(pred, target, loss_name="huber", huber_delta=1.0):

    if loss_name == "mse":
        return F.mse_loss(pred, target)

    elif loss_name == "mae":
        return F.l1_loss(pred, target)

    elif loss_name == "huber":
        return F.huber_loss(pred, target, delta=huber_delta)

    else:
        raise ValueError(f"Unknown loss: {loss_name}")

def compute_stage_loss(pred, target, stage_tasks):

    losses = []

    for t in stage_tasks:

        mask = ~torch.isnan(target[:, t])

        if mask.sum() == 0:
            continue

        pred_i = pred[mask, t]
        target_i = target[mask, t]

        loss_i = F.huber_loss(pred_i, target_i, delta=1.0)

        losses.append(loss_i)

    if len(losses) == 0:
        return torch.tensor(0.0, device=pred.device)

    return torch.stack(losses).mean()

########### PCGrad implmentation

class PCGrad:
    def __init__(self, optimizer):
        self.optimizer = optimizer

    @staticmethod
    def _flatten_grads(grads, params):
        flat_grads = []
        shapes = []

        for g, p in zip(grads, params):
            if g is None:
                g = torch.zeros_like(p)

            shapes.append(p.shape)
            flat_grads.append(g.reshape(-1))

        return torch.cat(flat_grads), shapes

    @staticmethod
    def _unflatten_grad(flat_grad, params):
        grads = []
        idx = 0

        for p in params:
            numel = p.numel()

            grads.append(
                flat_grad[idx:idx+numel].view_as(p)
            )

            idx += numel

        return grads

    def pc_backward(self, losses, params):

        valid_losses = [l for l in losses if l is not None]

        if len(valid_losses) == 0:
            return

        # ==========================================
        # Compute gradients with autograd.grad
        # ==========================================
        grad_vecs = []

        for i, loss in enumerate(valid_losses):

            retain = i < (len(valid_losses) - 1)

            grads = torch.autograd.grad(
                loss,
                params,
                retain_graph=retain,
                allow_unused=True,
            )

            flat_grad, _ = self._flatten_grads(grads, params)

            grad_vecs.append(flat_grad)

        G = torch.stack(grad_vecs)   # [T, P]

        # ==========================================
        # PCGrad projection
        # ==========================================
        proj_grads = G.clone()

        num_tasks = G.shape[0]

        for i in range(num_tasks):

            order = torch.randperm(num_tasks, device=G.device)

            for j in order:

                if i == j:
                    continue

                gij = torch.dot(proj_grads[i], G[j])

                if gij < 0:

                    gj_norm_sq = torch.dot(G[j], G[j]) + 1e-12

                    proj_grads[i] -= (
                        gij / gj_norm_sq
                    ) * G[j]

        # final_grad = proj_grads.mean(dim=0) # hurts training with pcgrad
        final_grad = proj_grads.sum(dim=0)

        # ==========================================
        # Write gradients back
        # ==========================================
        unflat = self._unflatten_grad(final_grad, params)

        self.optimizer.zero_grad(set_to_none=True)

        for p, g in zip(params, unflat):
            p.grad = g

    def step(self):

        torch.nn.utils.clip_grad_norm_(
            [
                p
                for group in self.optimizer.param_groups
                for p in group["params"]
            ],
            max_norm=1.0
        )

        self.optimizer.step()

# =========================
# TRAIN / EVAL EPOCHS
# =========================
def train_epoch(
    model,
    loader,
    optimizer,
    device,
    active_tasks=None,
    task_weights=None,
):

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
        out_detached = out.detach()
        # model._pcgrad_forward_cache = out   # cache forward for PCGrad
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

            max_tasks = getattr(model, "pcgrad_max_tasks", None)

            if max_tasks is not None and len(task_indices) > max_tasks:
                perm = torch.randperm(len(task_indices))
                task_indices = [task_indices[i] for i in perm[:max_tasks]]
            else:
                # still shuffle even if not truncating
                perm = torch.randperm(len(task_indices))
                task_indices = [task_indices[i] for i in perm]

            # ==========================================
            # Shared backbone params only
            # ==========================================

            shared_modules = [
                model.input_proj,
                model.convs,
                model.norms,
                model.shared,
                model.fp_mlp,
                model.global_proj,
            ]

            params = []

            for module in shared_modules:
                params.extend(
                    [p for p in module.parameters() if p.requires_grad]
                )

            task_losses_pcgrad = []

            for task_idx in task_indices:

                mask = ~torch.isnan(target[:, task_idx])

                if mask.sum() == 0:
                    continue

                pred_i = out[mask, task_idx]
                target_i = target[mask, task_idx]

                loss_i = compute_base_loss(
                    pred_i,
                    target_i,
                    loss_name=model.loss_name,
                    huber_delta=model.huber_delta
                )

                if task_weights is not None:
                    loss_i = loss_i * task_weights[task_idx]

                task_losses_pcgrad.append(loss_i)

            if len(task_losses_pcgrad) > 0:

                pcgrad.pc_backward(task_losses_pcgrad, params)

                torch.nn.utils.clip_grad_norm_(
                    model.parameters(),
                    max_norm=1.0
                )

                optimizer.step()

                # logging loss
                with torch.no_grad():
                    loss = torch.stack(
                        [l.detach() for l in task_losses_pcgrad]
                    ).mean()

            else:
                loss = torch.tensor(0.0, device=device)
        else:
            loss = multitask_loss(
                out,
                target,
                active_tasks=active_tasks,
                task_weights=task_weights,
            )

            if loss.requires_grad:
                optimizer.zero_grad(set_to_none=True)
                loss.backward()

                torch.nn.utils.clip_grad_norm_(
                    model.parameters(),
                    max_norm=1.0
                )

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

        all_preds.append(out_detached.cpu())
        all_labels.append(target.detach().cpu())

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
            out_detached = out.detach()
            target = batch.y.float()

            loss = multitask_loss(
                out,
                target,
                active_tasks=active_tasks
            )

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

def log_trainable_parameters(model, stage_idx):

    total = 0
    trainable = 0

    for p in model.parameters():
        total += p.numel()

        if p.requires_grad:
            trainable += p.numel()

    logger.info(
        f"Stage {stage_idx}: "
        f"{trainable:,}/{total:,} trainable "
        f"({100*trainable/total:.1f}%)"
    )

    for module_name in [
        "input_proj",
        "convs",
        "norms",
        "fp_mlp",
        "global_proj",
        "shared",
        "group_trunks",
    ]:

        module = getattr(model, module_name, None)

        if module is None:
            continue

        frozen = all(
            not p.requires_grad
            for p in module.parameters()
        )

        logger.debug(
            f"  {module_name}: "
            f"{'FROZEN' if frozen else 'TRAINABLE'}"
        )

def apply_curriculum_freezing(model, stage_idx, curriculum_cfg, context):

    strategy = curriculum_cfg.get(
        "strategy",
        "none"
    )

    # -------------------------
    # Default: train everything
    # -------------------------
    if strategy == "none":

        for p in model.parameters():
            p.requires_grad = True

        log_trainable_parameters(
            model,
            stage_idx
        )

        return
    
    elif strategy == "progressive_freezing_gradual":
        for p in model.parameters():
            p.requires_grad = True

        if stage_idx == 0:
            pass
        
        stage_name = curriculum_cfg["stages"][stage_idx]["name"]

        freeze_schedule = curriculum_cfg.get(
            "freeze_schedule",
            {}
        )

        if stage_name not in freeze_schedule:
            logger.warning(
                f"No freeze schedule found for stage {stage_name}"
            )

        modules_to_freeze = (
            freeze_schedule
            .get(stage_name, {})
            .get("freeze", [])
        )

        for module_name in modules_to_freeze:

            module = getattr(model, module_name, None)

            if module is None:
                logger.warning(
                    f"Unknown module in freeze schedule: {module_name}"
                )
                continue

            for p in module.parameters():
                p.requires_grad = False

        log_trainable_parameters(
            model,
            stage_idx
        )

        return
    
    elif strategy == "progressive_freezing_steep":

    # start with everything trainable
        for p in model.parameters():
            p.requires_grad = True

        if stage_idx >= 1:

            frozen_modules = [
                model.input_proj,
                model.convs,
                model.norms,
                model.shared,
                model.fp_mlp,
                model.global_proj,
            ]

            for module in frozen_modules:
                for p in module.parameters():
                    p.requires_grad = False
        
        log_trainable_parameters(
            model,
            stage_idx
        )

        return
    
    elif strategy == "representation_pretrain":

        frozen_backbone = sum(
            p.numel()
            for module in [
                model.input_proj,
                model.convs,
                model.norms,
                model.shared,
                model.fp_mlp,
                model.global_proj,
            ]
            for p in module.parameters()
            if not p.requires_grad
        )

        logger.info(
            f"Frozen backbone params: "
            f"{frozen_backbone:,}"
        )

        rep_end_name = curriculum_cfg.get(
            "representation_end_stage",
            "permeability"
        )

        stage_names = [
            s["name"]
            for s in curriculum_cfg["stages"]
        ]

        rep_end_idx = stage_names.index(
            rep_end_name
        )

        # everything trainable initially

        for p in model.parameters():
            p.requires_grad = True

        # after permeability:
        # freeze the learned representation

        if stage_idx > rep_end_idx:

            frozen_modules = [

                model.input_proj,
                model.convs,
                model.norms,
                model.shared,
                model.fp_mlp,
                model.global_proj,
            ]

            for module in frozen_modules:

                for p in module.parameters():

                    p.requires_grad = False

        log_trainable_parameters(
            model,
            stage_idx
        )

        return
    
    elif strategy == "progressive_unfreeze":
        # -------------------------
        # Freeze everything first
        # -------------------------
        for p in model.parameters():
            p.requires_grad = False

        # -------------------------
        # Identify modules
        # -------------------------
        backbone_blocks = [
            model.input_proj,
            model.convs,
            model.norms,
            model.shared,
            model.fp_mlp,
            model.global_proj,
        ]

        backbone_blocks = [
            m for m in backbone_blocks
            if m is not None
        ]

        # task heads = everything not in backbone
        backbone_param_ids = {
            id(p)
            for block in backbone_blocks
            for p in block.parameters()
        }

        head_params = [
            p
            for p in model.parameters()
            if id(p) not in backbone_param_ids
        ]

        # -------------------------
        # Always train heads
        # -------------------------
        for p in head_params:
            p.requires_grad = True

        # -------------------------
        # Progressive unfreezing
        # -------------------------
        if stage_idx >= 1:

            for p in model.shared.parameters():
                p.requires_grad = True

            for p in model.fp_mlp.parameters():
                p.requires_grad = True

            for p in model.global_proj.parameters():
                p.requires_grad = True

        if stage_idx >= 2:

            for p in model.convs.parameters():
                p.requires_grad = True

            for p in model.norms.parameters():
                p.requires_grad = True

        if stage_idx >= 3:

            for p in model.input_proj.parameters():
                p.requires_grad = True

        log_trainable_parameters(
            model,
            stage_idx
        )

        return
    
    else:
        raise ValueError(
            f"Unknown curriculum strategy: {strategy}"
        )

    # # -------------------------
    # # Logging
    # # -------------------------
    # trainable = sum(
    #     p.numel()
    #     for p in model.parameters()
    #     if p.requires_grad
    # )

    # total = sum(
    #     p.numel()
    #     for p in model.parameters()
    # )

    # logger.info(
    #     f"Stage {stage_idx}: "
    #     f"{trainable:,}/{total:,} params trainable "
    #     f"({100*trainable/total:.1f}%)"
    # )

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
    weighted_mse = curriculum_cfg.get("weighted_mse", False)
    reset_optimizer = curriculum_cfg.get("reset_optimizer", True)

    train_list = context["train_list"]
    val_list = context["val_list"]

   # Model init
    
    sample = data_list[0]

    group_architecture = {
        stage["name"]: {
            "hidden_dim": stage.get("hidden_dim", 256),
            "output_dim": stage.get("output_dim", 256),
        }
        for stage in curriculum_cfg["stages"]
    }

    logger.debug(
        f"Group architecture: {group_architecture}"
    )

    model = ModelClass(
        input_dim=sample.x.shape[1],
        edge_dim=sample.edge_attr.shape[1],
        global_feat_dim=sample.global_features.shape[-1],
        fp_dim=sample.fp.shape[-1],
        num_tasks=sample.y.shape[-1],
        task_groups=context.get("task_groups"),
        group_architecture=group_architecture,
        **{k: v for k, v in params.items() if k != "lr"}
    ).to(device)

    model.use_pcgrad = config.get("use_pcgrad", True)
    model.pcgrad_max_tasks = config.get("pcgrad_max_tasks", None)

    if model.use_pcgrad:
        logger.info("PCGrad ENABLED for training")
    else:
        logger.info("PCGrad DISABLED (using standard multitask loss)")

    best_loss = float("inf")
    best_model_state = None 

    logger.info(f"Curriculum freezing strategy applied: {strategy}")

    optimizer = None
    scheduler = None

    cumulative_tasks = []

    for stage_idx, stage in enumerate(stages):

        if isinstance(stage["tasks"][0], str):
            current_stage_tasks = [
                context["task_names"].index(t)
                for t in stage["tasks"]
            ]
        else:
            current_stage_tasks = stage["tasks"]

        # ADD NEW TASKS TO CUMULATIVE SET
        for t in current_stage_tasks:
            if t not in cumulative_tasks:
                cumulative_tasks.append(t)

        stage_train_list = filter_dataset_for_tasks(
            train_list,
            cumulative_tasks
        )

        stage_val_list = filter_dataset_for_tasks(
            val_list,
            cumulative_tasks
        )

        # =====================================
        # TASK WEIGHTS
        # =====================================
        task_weights = None

        if weighted_mse:
            task_weights = compute_task_weights(
                stage_train_list,
                num_tasks=sample.y.shape[-1],
            ).to(device)

            logger.info(
                f"Stage {stage['name']} task weights: "
                f"{task_weights.cpu().numpy().round(3)}"
            )

        # =====================================
        # REPLAY EARLIER TASKS
        # =====================================

        stage_val_list = filter_dataset_for_tasks(
            val_list,
            cumulative_tasks
        )

        logger.info(f"Stage {stage['name']} dataset size: train={len(stage_train_list)}, val={len(stage_val_list)}")

        if len(stage_train_list) == 0:
            raise RuntimeError(f"Empty train set for stage {stage['name']}")

        if len(stage_val_list) == 0:
            raise RuntimeError(f"Empty val set for stage {stage['name']}")

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

        patience = stage.get("patience", 10)

        logger.info(f"=== Curriculum stage: {stage['name']} | tasks={cumulative_tasks} ===")

        # History tracking
        train_losses = []
        val_losses = []
        train_task_hist = []
        val_task_hist = []

        apply_curriculum_freezing(model, stage_idx, curriculum_cfg, context)

        # =====================================
        # OPTIMISER RESET LOGIC
        # =====================================
        need_new_optimizer = (
            optimizer is None
            or reset_optimizer
        )

        if need_new_optimizer:

            backbone_modules = [
                model.input_proj,
                model.convs,
                model.norms,
                model.shared,
                model.fp_mlp,
                model.global_proj,
            ]

            backbone_params = []
            for module in backbone_modules:
                backbone_params.extend(
                [p for p in module.parameters()
                if p.requires_grad]
            )

            backbone_param_ids = set(id(p) for p in backbone_params)

            head_params = [
                p for p in model.parameters()
                if id(p) not in backbone_param_ids
                and p.requires_grad
            ]

            if strategy == "representation_pretrain":

                stage_backbone_lr = lr

            else:

                backbone_lr_decay = curriculum_cfg.get(
                    "backbone_lr_decay",
                    0.5
                )

                stage_backbone_lr = (
                    lr *
                    (backbone_lr_decay ** stage_idx)
                )

            optimizer = torch.optim.Adam(
                [
                    {
                        "params": backbone_params,
                        "lr": stage_backbone_lr,
                    },
                    {
                        "params": head_params,
                        "lr": lr,
                    }
                ]
            )

            # log the optimiser parameter count so we can verify optimiser is being reset correctly
 
            opt_params = sum(
                p.numel()
                for group in optimizer.param_groups
                for p in group["params"]
            )

            logger.debug(
                f"Optimizer contains {opt_params:,} parameters"
            )

            scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
                optimizer,
                mode="min",
                factor=0.5,
                patience=5,
                verbose=True
            )

            logger.info(
                f"Optimizer RESET at stage {stage['name']}"
            )

        else:
            # preserve optimizer state but refresh trainable params

            trainable_params = list(
                filter(lambda p: p.requires_grad, model.parameters())
            )

            existing_params = {
                p for group in optimizer.param_groups
                for p in group["params"]
            }

            new_params = [
                p for p in trainable_params
                if p not in existing_params
            ]

            if len(new_params) > 0:
                optimizer.add_param_group({
                    "params": new_params,
                    "lr": lr,
                })

                # Rebuild scheduler because param groups changed
                scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
                    optimizer,
                    mode="min",
                    factor=0.5,
                    patience=5,
                    verbose=True
                )

            logger.info(
                f"Optimizer PRESERVED at stage {stage['name']}"
            )
        
        stage_epochs = stage.get("epochs", config.get("max_epochs", 100))

        # Per-stage early stopping
        best_stage_loss = float("inf")
        best_stage_state = None
        counter = 0

        for epoch in range(stage_epochs):

            # =========================
            # TRAIN
            # =========================
            train_loss, train_preds, train_labels, train_rmse, train_r2 = train_epoch(
                model,
                train_loader,
                optimizer,
                device,
                active_tasks=cumulative_tasks,
                task_weights=task_weights,
            )

            # =========================
            # MASK TRAIN RMSE TO STAGE TASKS ONLY
            # =========================
            train_rmse_stage = train_rmse.copy()
            mask_arr = np.ones_like(train_rmse_stage, dtype=bool)
            mask_arr[cumulative_tasks] = False
            train_rmse_stage[mask_arr] = np.nan

            # =========================
            # VALIDATION
            # =========================

            val_loss, val_preds, val_labels, val_mse = eval_epoch(
                model,
                val_loader,
                device,
                active_tasks=cumulative_tasks
            )

            stage_val_loss = compute_stage_loss(
                val_preds.to(device),
                val_labels.to(device),
                cumulative_tasks
            ).item()

            scheduler.step(stage_val_loss)

            val_rmse = np.sqrt(val_mse)

            val_rmse_stage = val_rmse.copy()
            mask_arr = np.ones_like(val_rmse_stage, dtype=bool)
            mask_arr[cumulative_tasks] = False
            val_rmse_stage[mask_arr] = np.nan

            val_r2 = per_task_r2(
                val_preds.to(device),
                val_labels.to(device)
            ).detach().cpu().numpy()

            # =========================
            # STORE HISTORY
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
            if stage_val_loss < best_stage_loss:
                best_stage_loss = stage_val_loss
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
        val_loss, val_preds, val_labels, _ = eval_epoch(
            model,
            val_loader,
            device,
            active_tasks=cumulative_tasks
        )

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

        torch.save(
        {
            "model_state_dict": model.state_dict(),

            "input_dim": sample.x.shape[1],
            "edge_dim": sample.edge_attr.shape[1],
            "global_feat_dim": sample.global_features.shape[-1],
            "fp_dim": sample.fp.shape[-1],
            "num_tasks": sample.y.shape[-1],

            "params": params,

            "task_names": context["task_names"],

            "task_groups": context["task_groups"],

            "group_architecture":
                model.group_architecture,

            "label_scalers":
                context["label_scalers"],

            "label_transform_metadata":
                context["label_transform_metadata"],
        },
        model_path,
        )

        logger.info(f"Saved best model to {model_path}")

        # Store best model
        if stage_val_loss < best_loss:
            best_loss = stage_val_loss
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
    param_grid = config["param_grid"]
    training_mode = config.get("training_mode", "multitask")
    
    required_params = {
        "lr",
        "hidden_dim",
        "num_layers",
        "dropout",
    }

    missing = required_params - set(param_grid.keys())

    if missing:
        raise ValueError(
            "The grid search in the YAML is missing required hyperparameters: "
            f"{sorted(missing)}. You need to add a minimum of 1 value per hyperparameter to train the model."
            f"\n\nExample:"
            f"\n\nparam_grid:"
            f"\n  lr: [0.001]"
            f"\n  hidden_dim: [256]"
            f"\n  num_layers: [4]"
            f"\n  dropout: [0.1]"
        )

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

    train_list = context["train_list"]
    val_list = context["val_list"]

    train_loader = context["train_loader"]
    val_loader = context["val_loader"]

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

    task_names = context["task_names"]

    grouped_tasks = [
        t for group in config["task_groups"].values() for t in group
    ]

    missing = set(task_names) - set(grouped_tasks)
    extra = set(grouped_tasks) - set(task_names)

    if missing:
        raise ValueError(f"Tasks missing from grouping: {missing}")

    if extra:
        raise ValueError(f"Unknown tasks in grouping: {extra}")

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

            if train_loss is None:
                raise RuntimeError("train_loss is None")

            val_loss, val_preds, val_labels, val_mse = eval_epoch(model, val_loader, device)

            if val_loss is None:
                raise RuntimeError("val_loss is None")

            val_task_rmse = np.sqrt(val_mse)
            if val_loss is not None and not np.isnan(val_loss):
                scheduler.step(val_loss)

            train_losses.append(train_loss)
            val_losses.append(val_loss)

            train_task_hist.append(train_task_rmse)
            val_task_hist.append(val_task_rmse)

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

            scaler_path = config["scaler_path"]
            transform_metadata_path = config["transform_metadata_path"]

            Path(scaler_path).parent.mkdir(
                parents=True,
                exist_ok=True
            )

            joblib.dump(
                context["label_scalers"],
                scaler_path
            )

            with open(transform_metadata_path,"w") as f:
                json.dump(
                    context["label_transform_metadata"],
                    f,
                    indent=2
                )

            torch.save(
            {
                # -------------------------
                # model
                # -------------------------
                "model_state_dict": model.state_dict(),

                # architecture metadata
                "input_dim": input_dim,
                "edge_dim": edge_dim,
                "global_feat_dim": global_feat_dim,
                "fp_dim": fp_dim,
                "num_tasks": num_tasks,

                # hyperparameters
                "params": params,

                # task bookkeeping
                "task_names": task_names,
                "task_groups": context["task_groups"],

                # -------------------------
                # preprocessing (NEW)
                # -------------------------
                "label_scalers": context["label_scalers"],

                "label_transform_metadata":
                    context["label_transform_metadata"],

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