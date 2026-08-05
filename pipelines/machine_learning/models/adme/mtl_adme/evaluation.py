import numpy as np
# import torch
import matplotlib.pyplot as plt
import math
import pandas as pd
import seaborn as sns

from pathlib import Path
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.metrics import mean_absolute_error
from scipy.stats import spearmanr, pearsonr

from pipeline.logger import setup_logger
from models.adme.mtl_adme.ood import compute_ood_features

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

def load_split_smiles_from_context_or_disk(context, config):
    """
    Recover train and validation SMILES.

    Preferred:
      context["train_smiles"], context["val_smiles"]

    Fallback:
      predictions_dir/split_smiles.csv
    """

    if "train_smiles" in context and "val_smiles" in context:
        return context["train_smiles"], context["val_smiles"]

    pred_dir = Path(config.get("predictions_dir", "outputs/mtl_adme/preds"))
    split_path = pred_dir / "split_smiles.csv"

    if not split_path.exists():
        raise FileNotFoundError(
            f"Could not find split SMILES metadata at {split_path}. "
            "Make sure training.py calls save_split_metadata()."
        )

    split_df = pd.read_csv(split_path)

    train_smiles = (
        split_df
        .loc[split_df["split"] == "train", "smiles"]
        .astype(str)
        .tolist()
    )

    val_smiles = (
        split_df
        .loc[split_df["split"] == "val", "smiles"]
        .astype(str)
        .tolist()
    )

    return train_smiles, val_smiles


def build_prediction_long_table(
    y_true,
    y_pred,
    task_names,
    val_smiles,
    ood_df,
):
    """
    Long-format table:
      one row per molecule-task observation with a non-missing label.
    """

    rows = []

    for mol_idx, smi in enumerate(val_smiles):
        ood_row = ood_df.iloc[mol_idx].drop("smiles").to_dict()

        for task_idx, task in enumerate(task_names):
            yt = y_true[mol_idx, task_idx]

            if np.isnan(yt):
                continue

            yp = y_pred[mol_idx, task_idx]
            residual = yp - yt

            rows.append({
                "mol_index": mol_idx,
                "smiles": smi,
                "task": task,
                "y_true": yt,
                "y_pred": yp,
                "residual": residual,
                "abs_error": abs(residual),
                **ood_row,
            })

    return pd.DataFrame(rows)


def compute_metrics_by_ood_bin(pred_long_df):
    rows = []

    for task, task_df in pred_long_df.groupby("task"):
        for ood_bin, bin_df in task_df.groupby("ood_bin"):

            if len(bin_df) < 2:
                continue

            metrics = compute_metrics(
                bin_df["y_true"].values,
                bin_df["y_pred"].values,
            )

            rows.append({
                "task": task,
                "ood_bin": ood_bin,
                "n_samples": len(bin_df),
                **metrics,
            })

    return pd.DataFrame(rows)


def compute_overall_metrics_by_ood_bin(pred_long_df):
    rows = []

    for ood_bin, bin_df in pred_long_df.groupby("ood_bin"):

        if len(bin_df) < 2:
            continue

        metrics = compute_metrics(
            bin_df["y_true"].values,
            bin_df["y_pred"].values,
        )

        rows.append({
            "task": "__all_tasks__",
            "ood_bin": ood_bin,
            "n_samples": len(bin_df),
            **metrics,
        })

    return pd.DataFrame(rows)

def metrics_by_scaffold(pred_long_df):
    """
    Compare model performance on known vs novel Murcko scaffolds.
    """

    rows = []

    for known, df in pred_long_df.groupby("known_scaffold"):

        if len(df) < 5:
            continue

        metrics = compute_metrics(
            df["y_true"],
            df["y_pred"],
        )

        rows.append({
            "known_scaffold": known,
            "n_samples": len(df),
            **metrics
        })

    return pd.DataFrame(rows)

def metrics_by_cluster_size(
    pred_long_df,
    bins=(0,5,20,100,np.inf),
):
    """
    Evaluate prediction quality as a function of local
    chemical density.
    """

    df = pred_long_df.copy()

    labels = [
        "<5",
        "5-20",
        "20-100",
        ">100",
    ]

    df["cluster_bin"] = pd.cut(
        df["cluster_size"],
        bins=bins,
        labels=labels,
        include_lowest=True,
    )

    rows = []

    for cluster, sub in df.groupby("cluster_bin"):

        if len(sub) < 3:
            continue

        metrics = compute_metrics(
            sub["y_true"],
            sub["y_pred"],
        )

        rows.append({
            "cluster_bin": cluster,
            "n_samples": len(sub),
            **metrics
        })

    return pd.DataFrame(rows)

def metrics_by_similarity(
    pred_long_df,
):
    """
    Performance as a function of nearest-neighbour similarity.
    """

    df = pred_long_df.copy()

    bins = np.linspace(0,1,6)

    df["similarity_bin"] = pd.cut(
        df["nearest_train_tanimoto"],
        bins,
        include_lowest=True,
    )

    rows = []

    for sim, sub in df.groupby("similarity_bin"):

        if len(sub) < 5:
            continue

        metrics = compute_metrics(
            sub["y_true"],
            sub["y_pred"],
        )

        rows.append({
            "similarity_bin": str(sim),
            "n_samples": len(sub),
            **metrics
        })

    return pd.DataFrame(rows)

def metrics_by_physchem_distance(
    pred_long_df,
):
    """
    Performance versus physicochemical distance from the
    training distribution.
    """

    df = pred_long_df.copy()

    df["distance_bin"] = pd.qcut(

        df["physchem_robust_distance"],

        5,

        duplicates="drop",

    )

    rows = []

    for dist, sub in df.groupby("distance_bin"):

        if len(sub) < 5:
            continue

        metrics = compute_metrics(
            sub["y_true"],
            sub["y_pred"],
        )

        rows.append({

            "distance_bin": str(dist),

            "n_samples": len(sub),

            **metrics

        })

    return pd.DataFrame(rows)

def compute_error_correlations(pred_long_df):
    """
    Correlation between error and each OOD descriptor.
    """

    features = [
        "nearest_train_tanimoto",
        "mean_top5_train_tanimoto",
        "physchem_robust_distance",
        "cluster_size",
        "ood_score",
    ]

    rows = []

    for feature in features:

        df = pred_long_df[
            [feature, "abs_error"]
        ].dropna()

        if len(df) < 5:
            continue

        rho = spearmanr(
            df[feature],
            df["abs_error"],
        ).correlation

        r = pearsonr(
            df[feature],
            df["abs_error"],
        )[0]

        rows.append({

            "feature": feature,

            "pearson": r,

            "spearman": rho,

        })

    return pd.DataFrame(rows)

def plot_ood_diagnostics(pred_long_df, plot_dir):
    """
    Diagnostic plots:
      1. Absolute error vs nearest train similarity
      2. Absolute error vs OOD score
      3. MAE by OOD bin
    """

    MAX_POINTS = 20000

    if len(pred_long_df) > MAX_POINTS:
        plot_df = pred_long_df.sample(
            n=MAX_POINTS,
            random_state=42,
        )
    else:
        plot_df = pred_long_df

    ood_plot_dir = plot_dir / "ood"
    ood_plot_dir.mkdir(parents=True, exist_ok=True)

    # -------------------------
    # Overall abs error vs nearest similarity
    # -------------------------
    plt.figure(figsize=(7, 5))
    sns.scatterplot(
        data=plot_df,
        x="nearest_train_tanimoto",
        y="abs_error",
        hue="ood_bin",
        alpha=0.4,
        s=20,
    )
    
    sns.regplot(
        data=pred_long_df,
        x="nearest_train_tanimoto",
        y="abs_error",
        scatter=False,
        lowess=False,
        ci=False,
        color="black",
    )
    plt.xlabel("Nearest train Tanimoto")
    plt.ylabel("Absolute error")
    plt.title("Absolute error vs nearest train similarity")
    plt.tight_layout()
    plt.savefig(
        ood_plot_dir / "abs_error_vs_nearest_train_tanimoto.png",
        dpi=300
    )
    plt.close()

    # -------------------------
    # Overall abs error vs OOD score
    # -------------------------
    plt.figure(figsize=(7, 5))
    sns.scatterplot(
        data=plot_df,
        x="ood_score",
        y="abs_error",
        hue="ood_bin",
        alpha=0.4,
        s=20,
    )
    sns.regplot(
        data=pred_long_df,
        x="ood_score",
        y="abs_error",
        scatter=False,
        lowess=False,
        ci=False,
        color="black",
    )
    plt.xlabel("OOD score")
    plt.ylabel("Absolute error")
    plt.title("Absolute error vs OOD score")
    plt.tight_layout()
    plt.savefig(
        ood_plot_dir / "abs_error_vs_ood_score.png",
        dpi=300
    )
    plt.close()

    # -------------------------
    # MAE by OOD bin, all tasks pooled
    # -------------------------
    plt.figure(figsize=(6, 5))
    order = ["in_domain", "moderate_ood", "high_ood", "unknown"]

    sns.barplot(
        data=pred_long_df,
        x="ood_bin",
        y="abs_error",
        order=[x for x in order if x in pred_long_df["ood_bin"].unique()],
        estimator=np.mean,
        errorbar="se",
    )
    plt.xlabel("OOD bin")
    plt.ylabel("MAE")
    plt.title("MAE by OOD bin")
    plt.tight_layout()
    plt.savefig(
        ood_plot_dir / "mae_by_ood_bin.png",
        dpi=300
    )
    plt.close()


def plot_per_task_ood_panels(pred_long_df, plot_dir, task_names):
    ood_plot_dir = plot_dir / "ood"
    ood_plot_dir.mkdir(parents=True, exist_ok=True)

    tasks_present = [
        t for t in task_names
        if t in set(pred_long_df["task"])
    ]

    if len(tasks_present) == 0:
        return

    n_cols = 3
    n_rows = math.ceil(len(tasks_present) / n_cols)

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(5 * n_cols, 4 * n_rows)
    )

    axes = np.array(axes).flatten()

    for i, task in enumerate(tasks_present):
        ax = axes[i]
        task_df = pred_long_df[pred_long_df["task"] == task]

        if len(task_df) > 5000:
            task_plot_df = task_df.sample(
                n=5000,
                random_state=42,
            )
        else:
            task_plot_df = task_df

        if len(task_df) < 3:
            ax.axis("off")
            continue

        sns.scatterplot(
            data=task_plot_df,
            x="ood_score",
            y="abs_error",
            hue="ood_bin",
            alpha=0.5,
            s=18,
            ax=ax,
            legend=False,
        )

        try:
            sns.regplot(
                data=task_df,
                x="ood_score",
                y="abs_error",
                scatter=False,
                lowess=False,
                ci=False,
                color="black",
                ax=ax,
            )
        except Exception:
            pass

        ax.set_title(task)
        ax.set_xlabel("OOD score")
        ax.set_ylabel("Absolute error")
        ax.grid(True, alpha=0.3)

    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.savefig(
        ood_plot_dir / "per_task_abs_error_vs_ood_score.png",
        dpi=300
    )
    plt.close()

def plot_per_task_reliability_curves(
    task_reliability_df,
    plot_dir,
    max_cols=3,
):
    """
    Plot prediction reliability curves for every task.

    X-axis:
        OOD percentile (higher = less reliable)

    Y-axis:
        Mean absolute error

    Generates:
        1. Combined grid plot
        2. Individual task plots
    """

    reliability_dir = plot_dir / "ood" / "reliability_curves"
    reliability_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    if len(task_reliability_df) == 0:
        return


    # ======================================================
    # Individual task plots
    # ======================================================

    for task, task_df in task_reliability_df.groupby("task"):

        plt.figure(figsize=(6,5))

        task_df = task_df.sort_values(
            "mean_ood_percentile"
        )

        plt.plot(
            task_df["mean_ood_percentile"],
            task_df["mean_abs_error"],
            marker="o",
        )

        plt.xlabel(
            "OOD percentile\n(higher = less reliable)"
        )

        plt.ylabel(
            "Mean absolute error"
        )

        plt.title(
            task
        )

        plt.grid(True)

        plt.tight_layout()

        plt.savefig(
            reliability_dir / f"{task}_reliability_curve.png",
            dpi=300,
            bbox_inches="tight",
        )

        plt.close()


    # ======================================================
    # Combined grid
    # ======================================================

    tasks = task_reliability_df["task"].unique()

    n_tasks = len(tasks)

    n_rows = math.ceil(
        n_tasks / max_cols
    )

    fig, axes = plt.subplots(
        n_rows,
        max_cols,
        figsize=(
            5 * max_cols,
            4 * n_rows
        ),
    )

    axes = np.array(axes).flatten()


    for i, task in enumerate(tasks):

        ax = axes[i]

        task_df = (
            task_reliability_df[
                task_reliability_df["task"] == task
            ]
            .sort_values("mean_ood_percentile")
        )

        ax.plot(
            task_df["mean_ood_percentile"],
            task_df["mean_abs_error"],
            marker="o",
        )
        ax.set_title(task)
        ax.set_xlabel(
            "OOD percentile"
        )
        ax.set_ylabel(
            "MAE"
        )
        ax.grid(True)

    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])


    plt.tight_layout()

    plt.savefig(
        plot_dir /
        "ood" /
        "per_task_prediction_reliability_curves.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()

def plot_metric_table(
    df,
    category_col,
    plot_dir,
    filename,
    title,
):
    """
    Plot MAE and R² for grouped analyses.
    """

    plot_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    if len(df) == 0:
        return

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(10,4)
    )

    sns.barplot(
        data=df,
        x=category_col,
        y="mae",
        ax=axes[0],
    )

    axes[0].set_title("MAE")
    axes[0].tick_params(axis="x", rotation=30)

    sns.barplot(
        data=df,
        x=category_col,
        y="r2",
        ax=axes[1],
    )

    axes[1].set_title("R²")
    axes[1].tick_params(axis="x", rotation=30)

    fig.suptitle(title)

    plt.tight_layout()

    plt.savefig(
        plot_dir / filename,
        dpi=300,
    )

    plt.close()

def plot_error_correlations(
    corr_df,
    plot_dir,
):
    """
    Barplot of Pearson and Spearman correlation between
    OOD descriptors and prediction error.
    """

    if len(corr_df) == 0:
        return

    corr_df = corr_df.melt(
        id_vars="feature",
        value_vars=["pearson", "spearman"],
        var_name="metric",
        value_name="correlation",
    )

    plt.figure(figsize=(8,5))

    sns.barplot(
        data=corr_df,
        x="feature",
        y="correlation",
        hue="metric",
    )

    plt.xticks(rotation=30, ha="right")

    plt.tight_layout()

    plt.savefig(
        plot_dir / "ood_feature_correlations.png",
        dpi=300,
    )

    plt.close()

def compute_reliability_curve(
    pred_long_df,
    confidence_col="ood_percentile",
    n_bins=10,
):
    """
    Prediction reliability curve.

    Higher confidence percentile = more OOD.

    Reports:
        mean OOD percentile
        mean absolute error
        number of predictions
    """

    df = pred_long_df.copy()

    df["confidence_bin"] = pd.qcut(
        df[confidence_col],
        q=n_bins,
        duplicates="drop",
    )

    rows = []

    for b, sub in df.groupby("confidence_bin"):

        if len(sub) < 5:
            continue

        rows.append({
            "confidence_bin": str(b),
            "mean_ood_percentile":
                sub[confidence_col].mean(),
            "mean_abs_error":
                sub["abs_error"].mean(),
            "n_predictions":
                len(sub),
        })

    return pd.DataFrame(rows)

def plot_reliability_curve(
    reliability_df,
    plot_dir,
):

    plt.figure(figsize=(7,5))

    plt.plot(
        reliability_df["mean_ood_percentile"],
        reliability_df["mean_abs_error"],
        marker="o",
    )

    plt.xlabel(
        "OOD percentile (higher = less reliable)"
    )

    plt.ylabel(
        "Mean absolute error"
    )

    plt.title(
        "Prediction reliability curve"
    )

    plt.grid(True)

    plt.tight_layout()

    plt.savefig(
        plot_dir / "prediction_reliability_curve.png",
        dpi=300,
    )

    plt.close()

def compute_task_reliability_curves(
    pred_long_df,
    n_bins=10,
):

    rows = []

    for task, task_df in pred_long_df.groupby("task"):

        task_df = task_df.copy()

        task_df["confidence_bin"] = pd.qcut(
            task_df["ood_percentile"],
            q=n_bins,
            duplicates="drop",
        )

        for b, sub in task_df.groupby(
            "confidence_bin"
        ):

            if len(sub) < 5:
                continue

            rows.append({
                "task": task,
                "confidence_bin": str(b),
                "mean_ood_percentile":
                    sub["ood_percentile"].mean(),
                "mean_abs_error":
                    sub["abs_error"].mean(),
                "n_predictions":
                    len(sub),
            })

    return pd.DataFrame(rows)

def plot_task_reliability_curves(
    task_reliability_df,
    plot_dir,
):
    """
    Plot prediction reliability curves per task.
    """

    ood_plot_dir = plot_dir / "ood"
    ood_plot_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    plt.figure(figsize=(10, 7))

    for task, df in task_reliability_df.groupby("task"):

        plt.plot(
            df["mean_ood_percentile"],
            df["mean_abs_error"],
            marker="o",
            label=task,
            alpha=0.6,
        )

    plt.xlabel(
        "OOD percentile (higher = less reliable)"
    )

    plt.ylabel(
        "Mean absolute error"
    )

    plt.title(
        "Task-specific prediction reliability curves"
    )

    plt.legend(
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        fontsize=8,
    )

    plt.grid(True)

    plt.tight_layout()

    plt.savefig(
        ood_plot_dir /
        "task_prediction_reliability_curves.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()

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

        # metrics_by_task = {}

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

        # r2 = r2_score(y_t, y_p)
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

    # =========================
    # 6. OOD / APPLICABILITY-DOMAIN EVALUATION
    # =========================
    ood_cfg = config.get("ood", {})

    if ood_cfg.get("enabled", False):
        logger.info("Computing OOD/applicability-domain metrics")

        train_smiles, val_smiles = load_split_smiles_from_context_or_disk(
            context=context,
            config=config,
        )

        if len(val_smiles) != y_true.shape[0]:
            raise ValueError(
                "Number of validation SMILES does not match y_true rows: "
                f"len(val_smiles)={len(val_smiles)}, "
                f"y_true.shape[0]={y_true.shape[0]}. "
                "Check that validation predictions were generated with "
                "shuffle=False and that the SMILES metadata corresponds "
                "to this prediction run."
            )

        ood_df = compute_ood_features(
            query_smiles=val_smiles,
            train_smiles=train_smiles,
            config=ood_cfg,
        )

        ood_df.insert(
            0,
            "smiles",
            val_smiles
        )

        ood_df.to_csv(
            plot_dir / "ood_features_by_molecule.csv",
            index=False
        )

        pred_long_df = build_prediction_long_table(
            y_true=y_true,
            y_pred=y_pred,
            task_names=task_names,
            val_smiles=val_smiles,
            ood_df=ood_df,
        )

        # Convert raw OOD score into dataset-relative percentile
        # 0 = most in-domain
        # 100 = most OOD

        pred_long_df["ood_percentile"] = (
            pred_long_df["ood_score"]
            .rank(pct=True)
            * 100
        )

        pred_long_df.to_csv(
            plot_dir / "predictions_with_ood_scores_long.csv",
            index=False
        )

        metrics_by_ood_bin = compute_metrics_by_ood_bin(
            pred_long_df
        )

        overall_by_ood_bin = compute_overall_metrics_by_ood_bin(
            pred_long_df
        )

        metrics_by_ood_bin = pd.concat(
            [overall_by_ood_bin, metrics_by_ood_bin],
            ignore_index=True,
        )

        metrics_by_ood_bin.to_csv(
            plot_dir / "metrics_by_ood_bin.csv",
            index=False
        )

        scaffold_metrics = metrics_by_scaffold(pred_long_df)

        scaffold_metrics.to_csv(
            plot_dir / "metrics_by_scaffold.csv",
            index=False,
        )

        plot_metric_table(
            scaffold_metrics,
            "known_scaffold",
            plot_dir / "ood",
            "scaffold_metrics.png",
            "Known vs Novel Scaffolds",
        )

        cluster_metrics = metrics_by_cluster_size(pred_long_df)

        cluster_metrics.to_csv(
            plot_dir / "metrics_by_cluster_size.csv",
            index=False,
        )

        plot_metric_table(
            cluster_metrics,
            "cluster_bin",
            plot_dir / "ood",
            "cluster_density_metrics.png",
            "Prediction Quality vs Cluster Density",
        )

        similarity_metrics = metrics_by_similarity(pred_long_df)

        similarity_metrics.to_csv(
            plot_dir / "metrics_by_similarity.csv",
            index=False,
        )

        plot_metric_table(
            similarity_metrics,
            "similarity_bin",
            plot_dir / "ood",
            "similarity_metrics.png",
            "Prediction Quality vs Nearest Neighbour Similarity",
        )

        distance_metrics = metrics_by_physchem_distance(pred_long_df)

        distance_metrics.to_csv(
            plot_dir / "metrics_by_physchem_distance.csv",
            index=False,
        )

        plot_metric_table(

            distance_metrics,
            "distance_bin",
            plot_dir / "ood",
            "physchem_distance_metrics.png",
            "Prediction Quality vs Physicochemical Distance",
        )

        corr_df = compute_error_correlations(pred_long_df)

        reliability_df = compute_reliability_curve(
            pred_long_df
        )

        reliability_df.to_csv(
            plot_dir / "prediction_reliability_curve.csv",
            index=False,
        )

        task_reliability = compute_task_reliability_curves(
            pred_long_df
        )

        task_reliability.to_csv(
            plot_dir /
            "task_prediction_reliability_curves.csv",
            index=False,
        )

        plot_per_task_reliability_curves(
            task_reliability_df=task_reliability,
            plot_dir=plot_dir,
        )

        plot_task_reliability_curves(
            task_reliability,
            plot_dir,
        )

        plot_reliability_curve(
            reliability_df,
            plot_dir / "ood",
        )

        plot_error_correlations(
            corr_df,
            plot_dir / "ood",
        )

        reliability_df = compute_reliability_curve(
            pred_long_df,
            confidence_col="ood_percentile",
            n_bins=10,
        )

        reliability_df.to_csv(
            plot_dir / "prediction_reliability_curve.csv",
            index=False,
        )

        plot_reliability_curve(
            reliability_df,
            plot_dir / "ood",
        )

        corr_df.to_csv(
            plot_dir / "ood_feature_correlations.csv",
            index=False,
        )

        plot_ood_diagnostics(
            pred_long_df=pred_long_df,
            plot_dir=plot_dir,
        )

        plot_per_task_ood_panels(
            pred_long_df=pred_long_df,
            plot_dir=plot_dir,
            task_names=task_names,
        )

        logger.info(f"OOD evaluation written to {plot_dir/'ood'}")

    logger.info(f"Saved evaluation plots to {plot_dir}")

    return {}