import copy
import json

from pathlib import Path

import numpy as np
import torch
import matplotlib.pyplot as plt

from torch_geometric.loader import DataLoader

from pipeline.logger import setup_logger

from models.adme.mtl_adme.training import (
    custom_collate,
    train_epoch,
    eval_epoch,
    per_task_r2,
    compute_stage_loss,
    filter_dataset_for_tasks,
)

logger = setup_logger(__name__, debug_mode=False, simple_format=True)


def infer_active_tasks(data_list, task_names, min_label_count=5):
    num_tasks = len(task_names)
    counts = np.zeros(num_tasks, dtype=int)

    for data in data_list:
        y = data.y

        if y.dim() == 2:
            y = y.squeeze(0)

        y_np = y.detach().cpu().numpy()
        counts += np.isfinite(y_np).astype(int)

    active_tasks = [
        i for i, count in enumerate(counts)
        if count >= min_label_count
    ]

    logger.info(
        "Fine-tuning label counts: "
        + ", ".join(
            f"{task_names[i]}={counts[i]}"
            for i in range(num_tasks)
        )
    )

    logger.info(
        "Active fine-tuning tasks: "
        + ", ".join(task_names[i] for i in active_tasks)
    )

    return active_tasks, counts


def set_all_trainable(model, trainable):
    for p in model.parameters():
        p.requires_grad = trainable


def apply_finetune_freezing(model, freeze_mode):
    """
    Freeze modes
    ------------
    full:
        Train all parameters.

    last_layers:
        Freeze molecular encoders.
        Train shared trunk, group trunks, heads, gates/norms.

    heads_only:
        Freeze everything except task heads.

    group_heads:
        Freeze molecular encoders and shared trunk.
        Train group trunks and heads.
    """

    if freeze_mode == "full":
        set_all_trainable(model, True)

    elif freeze_mode == "last_layers":
        set_all_trainable(model, False)

        modules_to_train = [
            model.shared,
            model.group_trunks,
            model.heads,
            model.graph_norm,
            model.fp_norm,
            model.global_norm,
            model.pre_shared_norm,
        ]

        for module in modules_to_train:
            for p in module.parameters():
                p.requires_grad = True

        model.fp_gate.requires_grad = True
        model.global_gate.requires_grad = True

    elif freeze_mode == "heads_only":
        set_all_trainable(model, False)

        for p in model.heads.parameters():
            p.requires_grad = True

    elif freeze_mode == "group_heads":
        set_all_trainable(model, False)

        for p in model.group_trunks.parameters():
            p.requires_grad = True

        for p in model.heads.parameters():
            p.requires_grad = True

    else:
        raise ValueError(
            f"Unknown freeze_mode: {freeze_mode}. "
            "Use one of: full, last_layers, heads_only, group_heads."
        )

    total = sum(p.numel() for p in model.parameters())
    trainable = sum(
        p.numel()
        for p in model.parameters()
        if p.requires_grad
    )

    logger.info(
        f"Fine-tuning freeze_mode={freeze_mode}: "
        f"{trainable:,}/{total:,} trainable parameters "
        f"({100.0 * trainable / total:.2f}%)"
    )


def build_optimizer(model, lr, backbone_lr, weight_decay):
    backbone_module_names = [
        "input_proj",
        "convs",
        "norms",
        "fp_mlp",
        "global_proj",
        "film",
        "att_pool",
    ]

    backbone_params = []

    for name in backbone_module_names:
        module = getattr(model, name, None)

        if module is None:
            continue

        backbone_params.extend(
            [
                p for p in module.parameters()
                if p.requires_grad
            ]
        )

    backbone_param_ids = set(id(p) for p in backbone_params)

    head_params = [
        p for p in model.parameters()
        if p.requires_grad
        and id(p) not in backbone_param_ids
    ]

    param_groups = []

    if len(backbone_params) > 0:
        param_groups.append({
            "params": backbone_params,
            "lr": backbone_lr,
        })

    if len(head_params) > 0:
        param_groups.append({
            "params": head_params,
            "lr": lr,
        })

    if len(param_groups) == 0:
        raise RuntimeError(
            "No trainable parameters found. Check freeze_mode."
        )

    optimizer = torch.optim.AdamW(
        param_groups,
        weight_decay=weight_decay,
    )

    return optimizer


def save_finetuned_checkpoint(
    model,
    checkpoint,
    config,
    context,
    model_path,
    best_val_loss,
    best_params,
):
    model_path = Path(model_path)
    model_path.parent.mkdir(parents=True, exist_ok=True)

    sample = context["graphs"][0]

    torch.save(
        {
            "model_state_dict": model.state_dict(),

            "input_dim": sample.x.shape[1],
            "edge_dim": sample.edge_attr.shape[1],
            "global_feat_dim": sample.global_features.shape[-1],
            "fp_dim": sample.fp.shape[-1],
            "num_tasks": sample.y.shape[-1],

            "params": checkpoint["params"],

            "fine_tune_params": best_params,
            "fine_tuned_from": context.get("pretrained_checkpoint_path"),

            "task_names": checkpoint["task_names"],
            "task_groups": checkpoint["task_groups"],
            "group_architecture": checkpoint.get("group_architecture"),

            "label_scalers": checkpoint["label_scalers"],
            "label_transform_metadata": checkpoint["label_transform_metadata"],

            "global_feature_scaler": checkpoint.get("global_feature_scaler"),

            "best_val_loss": best_val_loss,
        },
        model_path,
    )

    logger.info(f"Saved fine-tuned model to {model_path}")


def fine_tune(context, config):
    device = context["device"]
    ModelClass = context["model_class"]

    checkpoint = context.get("adme_checkpoint")

    if checkpoint is None:
        checkpoint_path = config["checkpoint_path"]
        checkpoint = torch.load(
            checkpoint_path,
            map_location="cpu",
            weights_only=False,
        )
        context["adme_checkpoint"] = checkpoint
        context["pretrained_checkpoint_path"] = checkpoint_path

    ft_cfg = config.get("fine_tune", {})

    batch_size = ft_cfg.get("batch_size", config.get("batch_size", 16))
    max_epochs = ft_cfg.get("max_epochs", 50)
    patience = ft_cfg.get("patience", 10)
    lr = ft_cfg.get("lr", 1e-4)
    backbone_lr = ft_cfg.get("backbone_lr", lr * 0.1)
    weight_decay = ft_cfg.get("weight_decay", 1e-5)
    freeze_mode = ft_cfg.get("freeze_mode", "last_layers")
    min_label_count = ft_cfg.get("min_label_count_per_task", 5)

    train_list = context["train_list"]
    val_list = context["val_list"]

    task_names = checkpoint["task_names"]
    task_groups = checkpoint["task_groups"]

    if context.get("task_names") != task_names:
        raise ValueError(
            "Context task_names do not match checkpoint task_names.\n"
            f"Context: {context.get('task_names')}\n"
            f"Checkpoint: {task_names}"
        )

    sample = context["graphs"][0]

    model = ModelClass(
        input_dim=checkpoint["input_dim"],
        edge_dim=checkpoint["edge_dim"],
        global_feat_dim=checkpoint["global_feat_dim"],
        fp_dim=checkpoint["fp_dim"],
        num_tasks=checkpoint["num_tasks"],
        task_groups=task_groups,
        group_architecture=checkpoint.get("group_architecture"),
        **{
            k: v
            for k, v in checkpoint["params"].items()
            if k != "lr"
        },
    ).to(device)

    missing, unexpected = model.load_state_dict(
        checkpoint["model_state_dict"],
        strict=True,
    )

    if missing:
        raise RuntimeError(f"Missing checkpoint keys: {missing}")

    if unexpected:
        raise RuntimeError(f"Unexpected checkpoint keys: {unexpected}")

    logger.info("Loaded pretrained model weights for fine-tuning")

    model.loss_name = ft_cfg.get(
        "loss_name",
        getattr(model, "loss_name", "huber")
    )

    model.huber_delta = ft_cfg.get(
        "huber_delta",
        getattr(model, "huber_delta", 1.0)
    )

    model.use_pcgrad = False

    requested_active_tasks = ft_cfg.get("active_tasks", "auto")

    if requested_active_tasks == "auto":
        active_tasks, label_counts = infer_active_tasks(
            train_list,
            task_names,
            min_label_count=min_label_count,
        )
    else:
        active_tasks = [
            task_names.index(t)
            if isinstance(t, str)
            else int(t)
            for t in requested_active_tasks
        ]

        label_counts = None

    if len(active_tasks) == 0:
        raise RuntimeError(
            "No active fine-tuning tasks found. "
            "Check project label columns and min_label_count_per_task."
        )

    train_list = filter_dataset_for_tasks(
        train_list,
        active_tasks,
    )

    val_list = filter_dataset_for_tasks(
        val_list,
        active_tasks,
    )

    if len(train_list) == 0:
        raise RuntimeError(
            "Fine-tuning train set is empty after filtering for active tasks."
        )

    if len(val_list) == 0:
        raise RuntimeError(
            "Fine-tuning validation set is empty after filtering for active tasks."
        )

    train_loader = DataLoader(
        train_list,
        batch_size=batch_size,
        shuffle=True,
        collate_fn=custom_collate,
    )

    val_loader = DataLoader(
        val_list,
        batch_size=batch_size,
        shuffle=False,
        collate_fn=custom_collate,
    )

    apply_finetune_freezing(
        model,
        freeze_mode=freeze_mode,
    )

    optimizer = build_optimizer(
        model,
        lr=lr,
        backbone_lr=backbone_lr,
        weight_decay=weight_decay,
    )

    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode="min",
        factor=0.5,
        patience=5,
        verbose=True,
    )

    loss_curve_dir = Path(
        config.get(
            "loss_curve_dir",
            "outputs/mtl_adme_finetune/loss_curves"
        )
    )
    loss_curve_dir.mkdir(parents=True, exist_ok=True)

    pred_dir = Path(
        config.get(
            "predictions_dir",
            "outputs/mtl_adme_finetune/preds"
        )
    )
    pred_dir.mkdir(parents=True, exist_ok=True)

    model_path = Path(
        config.get(
            "model_path",
            "outputs/mtl_adme_finetune/models/mtl_adme_project_ft_best.pth"
        )
    )

    best_val_loss = float("inf")
    best_state = copy.deepcopy(
        {
            k: v.detach().cpu()
            for k, v in model.state_dict().items()
        }
    )

    counter = 0

    train_losses = []
    val_losses = []
    train_task_hist = []
    val_task_hist = []

    for epoch in range(max_epochs):
        train_loss, train_preds, train_labels, train_rmse, train_r2 = train_epoch(
            model,
            train_loader,
            optimizer,
            device,
            active_tasks=active_tasks,
            task_weights=None,
        )

        val_loss, val_preds, val_labels, val_mse = eval_epoch(
            model,
            val_loader,
            device,
            active_tasks=active_tasks,
        )

        stage_val_loss = compute_stage_loss(
            val_preds.to(device),
            val_labels.to(device),
            active_tasks,
        ).item()

        scheduler.step(stage_val_loss)

        val_rmse = np.sqrt(val_mse)

        train_rmse_stage = train_rmse.copy()
        val_rmse_stage = val_rmse.copy()

        inactive_mask = np.ones_like(train_rmse_stage, dtype=bool)
        inactive_mask[active_tasks] = False

        train_rmse_stage[inactive_mask] = np.nan
        val_rmse_stage[inactive_mask] = np.nan

        train_losses.append(train_loss)
        val_losses.append(val_loss)
        train_task_hist.append(train_rmse_stage)
        val_task_hist.append(val_rmse_stage)

        val_r2 = per_task_r2(
            val_preds.to(device),
            val_labels.to(device),
        ).detach().cpu().numpy()

        logger.info(
            f"Fine-tune epoch {epoch} | "
            f"Train Loss {train_loss:.4f} | "
            f"Val Loss {val_loss:.4f} | "
            f"Active Val Loss {stage_val_loss:.4f}"
        )

        logger.debug(f"Fine-tune Val RMSE: {np.round(val_rmse_stage, 3)}")
        logger.debug(f"Fine-tune Val R2:   {np.round(val_r2, 3)}")

        if stage_val_loss < best_val_loss:
            best_val_loss = stage_val_loss
            best_state = copy.deepcopy(
                {
                    k: v.detach().cpu()
                    for k, v in model.state_dict().items()
                }
            )
            counter = 0
        else:
            counter += 1

        if counter >= patience:
            logger.info("Fine-tuning early stopping triggered")
            break

    model.load_state_dict(best_state)
    model.to(device)

    final_val_loss, val_preds, val_labels, val_mse = eval_epoch(
        model,
        val_loader,
        device,
        active_tasks=active_tasks,
    )

    y_pred = val_preds.cpu().numpy()
    y_true = val_labels.cpu().numpy()

    np.save(pred_dir / "y_pred.npy", y_pred)
    np.save(pred_dir / "y_true.npy", y_true)

    context["y_pred"] = y_pred
    context["y_true"] = y_true

    train_task_hist = np.array(train_task_hist)
    val_task_hist = np.array(val_task_hist)

    plt.figure()
    plt.plot(train_losses, label="Train")
    plt.plot(val_losses, label="Val")
    plt.legend()
    plt.title("Fine-tuning loss")
    plt.savefig(loss_curve_dir / "fine_tuning_overall_loss.png", dpi=150)
    plt.close()

    n_tasks = len(task_names)
    n_cols = 3
    n_rows = int(np.ceil(n_tasks / n_cols))

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(5 * n_cols, 4 * n_rows),
    )

    axes = axes.flatten()

    for i in range(n_tasks):
        axes[i].plot(train_task_hist[:, i], label="Train")
        axes[i].plot(val_task_hist[:, i], label="Val")
        axes[i].set_title(task_names[i])
        axes[i].legend()

    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.savefig(loss_curve_dir / "fine_tuning_per_task_rmse.png", dpi=150)
    plt.close()

    best_params = {
        "freeze_mode": freeze_mode,
        "lr": lr,
        "backbone_lr": backbone_lr,
        "weight_decay": weight_decay,
        "max_epochs": max_epochs,
        "patience": patience,
        "active_tasks": [
            task_names[i]
            for i in active_tasks
        ],
    }

    save_finetuned_checkpoint(
        model=model,
        checkpoint=checkpoint,
        config=config,
        context=context,
        model_path=model_path,
        best_val_loss=best_val_loss,
        best_params=best_params,
    )

    return {
        "model": model,
        "best_val_loss": best_val_loss,
        "best_params": best_params,
        "task_names": task_names,
        "active_tasks": active_tasks,
        "predictions_dir": str(pred_dir),
        "model_path": str(model_path),
        "y_pred": y_pred,
        "y_true": y_true,
    }