import numpy as np
import torch
import matplotlib.pyplot as plt
import math
import pandas as pd
import seaborn as sns

from pathlib import Path
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.metrics import mean_absolute_error
from scipy.stats import spearmanr, pearsonr

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

#### HELPERS

def compute_metrics(y_true, y_pred):

    r2 = r2_score(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    mae = mean_absolute_error(y_true, y_pred)
    median_ae = np.median(np.abs(y_true - y_pred))

    try:
        spearman = spearmanr(y_true, y_pred).correlation
    except Exception:
        spearman = np.nan

    try:
        pearson = pearsonr(y_true, y_pred)[0]
    except Exception:
        pearson = np.nan

    return {
        "r2": r2,
        "rmse": rmse,
        "mae": mae,
        "median_ae": median_ae,
        "spearman": spearman,
        "pearson": pearson,
    }

def evaluate(context, config):
    plot_dir = Path(config.get("eval_plot_dir", "outputs/eval"))
    plot_dir.mkdir(parents=True, exist_ok=True)

    # -------------------------
    # LOAD DATA
    # -------------------------
    if "y_pred" in context and "y_true" in context:
        y_pred = context["y_pred"]
        y_true = context["y_true"]

        logger.info("Using predictions from context")

    else:
        pred_dir = Path(config.get("predictions_dir", "outputs/mtl_adme/preds"))

        y_pred = np.load(pred_dir / "y_pred.npy")
        y_true = np.load(pred_dir / "y_true.npy")

        logger.info(f"Loaded predictions from {pred_dir}")

    task_names = context.get("task_names", None)

    # --- Fix shapes ---
    if y_true.ndim == 3:
        y_true = np.squeeze(y_true, axis=1)
    if y_pred.ndim == 3:
        y_pred = np.squeeze(y_pred, axis=1)
    
    # =========================
    # INVERSE TRANSFORM
    # =========================
    scalers = context.get("label_scalers", None)

    if scalers is not None:
        logger.info("Applying inverse label scaling")

        metrics_by_task = {}

        for i, task in enumerate(task_names):
            scaler = scalers.get(task)

            if scaler is None:
                continue

            mask = ~np.isnan(y_true[:, i])

            if mask.sum() == 0:
                continue

            y_true[mask, i] = scaler.inverse_transform(
                y_true[mask, i].reshape(-1, 1)
            ).flatten()

            y_pred[mask, i] = scaler.inverse_transform(
                y_pred[mask, i].reshape(-1, 1)
            ).flatten()

    n_tasks = y_true.shape[1]

    if task_names is None:
        task_names = [f"Task_{i}" for i in range(n_tasks)]

    # =========================
    # 1. METRICS + COUNTS CSV
    # =========================
    rows = []

    metrics_by_task = {}

    for i, task in enumerate(task_names):

        mask = ~np.isnan(y_true[:, i])

        if mask.sum() == 0:
            continue

        y_t = y_true[mask, i]
        y_p = y_pred[mask, i]

        metrics = compute_metrics(y_t, y_p)

        metrics_by_task[task] = metrics

        logger.info(
            f"{task}: "
            f"R2={metrics['r2']:.3f}, "
            f"RMSE={metrics['rmse']:.3f}, "
            f"MAE={metrics['mae']:.3f}"
        )

        rows.append({
            "task": task,
            "n_samples": int(mask.sum()),
            **metrics
        })

    pd.DataFrame(rows).to_csv(plot_dir / "task_metrics.csv", index=False)

    # =========================
    # 2. TASK CORRELATION HEATMAP
    # =========================
    df_tasks = pd.DataFrame(y_true, columns=task_names)

    plt.figure(figsize=(6, 5))
    sns.heatmap(df_tasks.corr(), annot=True, cmap="coolwarm_r", center=0)
    plt.title("Task correlation matrix")
    plt.tight_layout()
    plt.savefig(plot_dir / "task_correlation.png", dpi=300)
    plt.close()

    # =========================
    # 3. TRUE VS PREDICTED
    # =========================
    n_cols = 3
    n_rows = math.ceil(n_tasks / n_cols)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 5*n_rows))
    axes = axes.flatten()

    for i, task in enumerate(task_names):
        ax = axes[i]

        mask = ~np.isnan(y_true[:, i])
        if mask.sum() == 0:
            ax.set_title(f"{task} (no data)")
            ax.axis('off')
            continue

        y_t = y_true[mask, i]
        y_p = y_pred[mask, i]

        ax.scatter(y_t, y_p, alpha=0.5)

        min_val = min(y_t.min(), y_p.min())
        max_val = max(y_t.max(), y_p.max())
        ax.plot([min_val, max_val], [min_val, max_val], 'r--')

        r2 = r2_score(y_t, y_p)
        metrics = metrics_by_task[task]

        ax.set_title(
            f"{task}\n"
            f"R²={metrics['r2']:.2f} | "
            f"RMSE={metrics['rmse']:.2f}\n"
            f"MAE={metrics['mae']:.2f} | "
            f"ρ={metrics['spearman']:.2f}"
        )
        ax.set_xlabel("True")
        ax.set_ylabel("Predicted")
        ax.grid(True)

    for j in range(i+1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.savefig(plot_dir / "true_vs_pred.png", dpi=300)
    plt.close()

    # =========================
    # 4. RESIDUALS
    # =========================
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 5*n_rows))
    axes = axes.flatten()

    for i, task in enumerate(task_names):
        ax = axes[i]

        mask = ~np.isnan(y_true[:, i])
        if mask.sum() == 0:
            ax.set_title(f"{task} (no data)")
            ax.axis('off')
            continue

        y_t = y_true[mask, i]
        y_p = y_pred[mask, i]

        residuals = y_p - y_t

        ax.scatter(y_t, residuals, alpha=0.5)
        ax.axhline(0, linestyle='--')

        ax.set_title(task)
        ax.set_xlabel("True")
        ax.set_ylabel("Residual")
        ax.grid(True)

    for j in range(i+1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.savefig(plot_dir / "residuals.png", dpi=300)
    plt.close()

    # =========================
    # 5. RESIDUAL HISTOGRAMS
    # =========================

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    axes = axes.flatten()

    for i, task in enumerate(task_names):

        ax = axes[i]

        mask = ~np.isnan(y_true[:, i])

        if mask.sum() == 0:
            ax.axis("off")
            continue

        y_t = y_true[mask, i]
        y_p = y_pred[mask, i]

        residuals = y_p - y_t

        ax.hist(residuals, bins=40)
        ax.set_title(task)
        ax.set_xlabel("Residual")
        ax.set_ylabel("Count")

    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.savefig(plot_dir / "residual_histograms.png", dpi=300)
    plt.close()

    logger.info(f"Saved evaluation plots to {plot_dir}")

    return {}