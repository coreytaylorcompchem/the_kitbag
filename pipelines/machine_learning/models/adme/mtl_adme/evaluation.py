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

palette = {
    "in_domain": "#7C3AED",
    "moderate_ood": "#1D4ED8",
    "high_ood": "#0F766E",
    "extreme_ood": "#92400E",
    "unknown": "#9CA3AF",
}

hue_order = [
    "in_domain",
    "moderate_ood",
    "high_ood",
    "extreme_ood",
    "unknown",
]

#### HELPERS

def metrics_by_molecule_group(
    y_true,
    y_pred,
    task_names,
    group_values,
    group_name,
    min_n=5,
):
    """
    Compute per-task metrics for molecule-level groups without building
    a molecule-task long dataframe.

    Parameters
    ----------
    group_values : array-like, shape [n_molecules]
        One group label per molecule, e.g. ood_bin or known_scaffold.
    """

    rows = []

    group_values = np.asarray(group_values)

    for task_idx, task in enumerate(task_names):
        task_true = y_true[:, task_idx]
        task_pred = y_pred[:, task_idx]

        valid_task = ~np.isnan(task_true)

        if valid_task.sum() == 0:
            continue

        for group in pd.unique(group_values):
            group_mask = group_values == group
            mask = valid_task & group_mask

            if mask.sum() < min_n:
                continue

            metrics = compute_metrics(
                task_true[mask],
                task_pred[mask],
            )

            rows.append({
                "task": task,
                group_name: group,
                "n_samples": int(mask.sum()),
                **metrics,
            })

    return pd.DataFrame(rows)

def overall_metrics_by_molecule_group(
    y_true,
    y_pred,
    group_values,
    group_name,
    min_n=20,
):
    """
    Compute pooled all-task metrics by molecule-level group without
    constructing a long dataframe.

    This loops through groups and flattens only the relevant block.
    """

    rows = []

    group_values = np.asarray(group_values)

    for group in pd.unique(group_values):
        mol_mask = group_values == group

        if mol_mask.sum() == 0:
            continue

        yt = y_true[mol_mask, :]
        yp = y_pred[mol_mask, :]

        valid = ~np.isnan(yt)

        if valid.sum() < min_n:
            continue

        metrics = compute_metrics(
            yt[valid],
            yp[valid],
        )

        rows.append({
            "task": "__all_tasks__",
            group_name: group,
            "n_samples": int(valid.sum()),
            **metrics,
        })

    return pd.DataFrame(rows)

def metrics_by_continuous_molecule_feature(
    y_true,
    y_pred,
    task_names,
    values,
    feature_name,
    bins,
    labels=None,
    min_n=5,
):
    """
    Per-task metrics by bins of a molecule-level continuous feature.
    """

    values = np.asarray(values, dtype=float)

    if labels is None:
        labels = [str(i) for i in range(len(bins) - 1)]

    binned = pd.cut(
        values,
        bins=bins,
        labels=labels,
        include_lowest=True,
    )

    per_task = metrics_by_molecule_group(
        y_true=y_true,
        y_pred=y_pred,
        task_names=task_names,
        group_values=np.asarray(binned).astype(object),
        group_name=f"{feature_name}_bin",
        min_n=min_n,
    )

    overall = overall_metrics_by_molecule_group(
        y_true=y_true,
        y_pred=y_pred,
        group_values=np.asarray(binned).astype(object),
        group_name=f"{feature_name}_bin",
        min_n=min_n,
    )

    return pd.concat(
        [overall, per_task],
        ignore_index=True,
    )

def metrics_by_quantile_molecule_feature(
    y_true,
    y_pred,
    task_names,
    values,
    feature_name,
    q=5,
    min_n=5,
):
    """
    Metrics by quantile bins of a molecule-level continuous feature.
    """

    values = pd.Series(values)

    binned = pd.qcut(
        values,
        q=q,
        duplicates="drop",
    )

    per_task = metrics_by_molecule_group(
        y_true=y_true,
        y_pred=y_pred,
        task_names=task_names,
        group_values=np.asarray(binned).astype(object),
        group_name=f"{feature_name}_bin",
        min_n=min_n,
    )

    overall = overall_metrics_by_molecule_group(
        y_true=y_true,
        y_pred=y_pred,
        group_values=np.asarray(binned).astype(object),
        group_name=f"{feature_name}_bin",
        min_n=min_n,
    )

    return pd.concat(
        [overall, per_task],
        ignore_index=True,
    )

def build_prediction_sample_table(
    y_true,
    y_pred,
    task_names,
    val_smiles,
    ood_df,
    task_error_scales=None,
    max_rows=100000,
    random_seed=42,
):
    """
    Build a sampled long-format prediction table for plotting only.

    This avoids materialising all molecule-task observations.
    """

    rng = np.random.default_rng(random_seed)

    valid_pairs = np.argwhere(~np.isnan(y_true))

    if len(valid_pairs) == 0:
        return pd.DataFrame()

    if len(valid_pairs) > max_rows:
        chosen = rng.choice(
            len(valid_pairs),
            size=max_rows,
            replace=False,
        )
        valid_pairs = valid_pairs[chosen]

    rows = []

    for mol_idx, task_idx in valid_pairs:
        yt = y_true[mol_idx, task_idx]
        yp = y_pred[mol_idx, task_idx]
        residual = yp - yt

        ood_row = ood_df.iloc[int(mol_idx)].to_dict()

        task = task_names[int(task_idx)]
        abs_error = float(abs(residual))

        if task_error_scales is not None:
            scale = task_error_scales.get(task, np.nan)

            if np.isfinite(scale) and scale > 0:
                normalized_abs_error = abs_error / scale
            else:
                normalized_abs_error = np.nan
        else:
            normalized_abs_error = np.nan

        rows.append({
            "mol_index": int(mol_idx),
            "smiles": val_smiles[int(mol_idx)],
            "task": task,
            "y_true": float(yt),
            "y_pred": float(yp),
            "residual": float(residual),
            "abs_error": abs_error,
            "normalized_abs_error": normalized_abs_error,
            **ood_row,
        })

    return pd.DataFrame(rows)

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

def compute_task_error_scales(
    y_true,
    y_pred,
    task_names,
    method="mae",
    min_scale=1e-8,
):
    """
    Compute per-task error scaling factors.

    Used to normalize absolute errors across endpoints.

    method="mae":
        scale = global validation MAE for that task

    method="iqr":
        scale = IQR of y_true for that task
    """

    scales = {}

    for task_idx, task in enumerate(task_names):
        yt = y_true[:, task_idx]
        yp = y_pred[:, task_idx]

        mask = ~np.isnan(yt)

        if mask.sum() < 2:
            scales[task] = np.nan
            continue

        if method == "mae":
            scale = np.mean(
                np.abs(yp[mask] - yt[mask])
            )

        elif method == "iqr":
            q75, q25 = np.nanpercentile(
                yt[mask],
                [75, 25],
            )
            scale = q75 - q25

        else:
            raise ValueError(
                f"Unknown task error scale method: {method}"
            )

        if not np.isfinite(scale) or scale < min_scale:
            scale = np.nan

        scales[task] = scale

    return scales

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

def compute_error_correlations(
    pred_long_df,
    candidate_features=None,
):
    """
    Correlation between absolute error and available OOD descriptors.

    This version is robust to optional OOD features being absent, e.g.
    cluster_size after removing Butina clustering for large datasets.
    """

    if candidate_features is None:
        candidate_features = [
            "nearest_train_tanimoto",
            "mean_top5_train_tanimoto",
            "physchem_robust_distance",
            "n_train_neighbors_ge_0.4",
            "n_train_neighbors_ge_0.5",
            "n_train_neighbors_ge_0.6",
            "n_train_neighbors_ge_0.7",
            "fraction_train_neighbors_ge_0.4",
            "fraction_train_neighbors_ge_0.5",
            "fraction_train_neighbors_ge_0.6",
            "fraction_train_neighbors_ge_0.7",
            "ood_score",
            "ood_percentile",
        ]

    rows = []

    available_features = [
        f for f in candidate_features
        if f in pred_long_df.columns
    ]

    missing_features = [
        f for f in candidate_features
        if f not in pred_long_df.columns
    ]

    if missing_features:
        logger.debug(
            "Skipping unavailable OOD correlation features: "
            f"{missing_features}"
        )

    for feature in available_features:
        df = pred_long_df[
            [feature, "abs_error"]
        ].dropna()

        if len(df) < 5:
            continue

        try:
            rho = spearmanr(
                df[feature],
                df["abs_error"],
            ).correlation
        except Exception:
            rho = np.nan

        try:
            r = pearsonr(
                df[feature],
                df["abs_error"],
            )[0]
        except Exception:
            r = np.nan

        rows.append({
            "feature": feature,
            "pearson": r,
            "spearman": rho,
            "n_samples": len(df),
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
        hue_order=hue_order,
        palette=palette,
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
        hue_order=hue_order,
        palette=palette,
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
    # Normalized abs error vs OOD score
    # -------------------------
    if "normalized_abs_error" in pred_long_df.columns:
        plt.figure(figsize=(7, 5))

        sns.scatterplot(
            data=plot_df,
            x="ood_score",
            y="normalized_abs_error",
            hue="ood_bin",
            hue_order=hue_order,
            palette=palette,
            alpha=0.4,
            s=20,
        )

        sns.regplot(
            data=pred_long_df.dropna(
                subset=["ood_score", "normalized_abs_error"]
            ),
            x="ood_score",
            y="normalized_abs_error",
            scatter=False,
            lowess=False,
            ci=False,
            color="black",
        )

        plt.xlabel("OOD score")
        plt.ylabel("Normalized absolute error")
        plt.title("Normalized absolute error vs OOD score")
        plt.tight_layout()

        plt.savefig(
            ood_plot_dir / "normalized_abs_error_vs_ood_score.png",
            dpi=300,
        )

        plt.close()
    
    # -------------------------
    # Normalized MAE by OOD bin
    # -------------------------
    if "normalized_abs_error" in pred_long_df.columns:
        plt.figure(figsize=(6, 5))

        order = [
            "in_domain",
            "moderate_ood",
            "high_ood",
            "extreme_ood",
            "unknown",
        ]

        sns.barplot(
            data=pred_long_df,
            x="ood_bin",
            y="normalized_abs_error",
            order=[
                x for x in order
                if x in pred_long_df["ood_bin"].unique()
            ],
            estimator=np.mean,
            errorbar="se",
            palette=palette,
        )

        plt.xlabel("OOD bin")
        plt.ylabel("Normalized MAE")
        plt.title("Task-normalized MAE by OOD bin")
        plt.tight_layout()

        plt.savefig(
            ood_plot_dir / "normalized_mae_by_ood_bin.png",
            dpi=300,
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
            hue_order=hue_order,
            palette=palette,            
            alpha=0.5,
            s=18,
            ax=ax,
            legend=(i == 0),
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

        # handles, labels = axes[0].get_legend_handles_labels()

        # fig.legend(
        #     handles,
        #     labels,
        #     loc="lower center",
        #     bbox_to_anchor=(0.5, -0.02),
        #     ncol=len(labels),
        #     frameon=False,
        # )
        
        ax.set_title(task)
        ax.set_xlabel("OOD score")
        ax.set_ylabel("Absolute error")
        ax.grid(True, alpha=0.3)

    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])
    
    handles, labels = axes[0].get_legend_handles_labels()

    if len(labels) > 0:
        fig.legend(
            handles,
            labels,
            loc="lower center",
            bbox_to_anchor=(0.5, 0.01),
            ncol=len(labels),
            frameon=False,
            title="OOD bin",
        )

    plt.tight_layout(rect=[0, 0.06, 1, 1])

    plt.savefig(
        ood_plot_dir / "per_task_abs_error_vs_ood_score.png",
        dpi=300,
        bbox_inches="tight",
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
    available OOD descriptors and prediction error.
    """

    if corr_df is None or len(corr_df) == 0:
        logger.info("No OOD feature correlations to plot")
        return

    required_cols = {
        "feature",
        "pearson",
        "spearman",
    }

    missing = required_cols - set(corr_df.columns)

    if missing:
        logger.warning(
            "Cannot plot OOD feature correlations because columns are missing: "
            f"{sorted(missing)}"
        )
        return

    corr_long = corr_df.melt(
        id_vars="feature",
        value_vars=["pearson", "spearman"],
        var_name="metric",
        value_name="correlation",
    )

    plt.figure(figsize=(8, 5))

    sns.barplot(
        data=corr_long,
        x="feature",
        y="correlation",
        hue="metric",
    )

    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()

    plot_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    plt.savefig(
        plot_dir / "ood_feature_correlations.png",
        dpi=300,
    )

    plt.close()

def compute_reliability_curve(
    pred_long_df,
    confidence_col="ood_percentile",
    error_col="abs_error",
    n_bins=10,
):
    """
    Prediction reliability curve.

    X-axis:
        OOD percentile or other confidence/risk column.

    Y-axis:
        Mean error column, e.g. abs_error or normalized_abs_error.
    """

    df = pred_long_df[
        [confidence_col, error_col]
    ].dropna().copy()

    if len(df) < 5:
        return pd.DataFrame()

    df["confidence_bin"] = pd.qcut(
        df[confidence_col],
        q=n_bins,
        duplicates="drop",
    )

    rows = []

    for b, sub in df.groupby("confidence_bin", observed=False):

        if len(sub) < 5:
            continue

        rows.append({
            "confidence_bin": str(b),
            "mean_ood_percentile": sub[confidence_col].mean(),
            f"mean_{error_col}": sub[error_col].mean(),
            "n_predictions": len(sub),
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

def plot_reliability_curve_generic(
    reliability_df,
    plot_dir,
    y_col,
    filename,
    ylabel,
    title,
):
    """
    Plot a reliability curve using a specified y column.
    """

    if reliability_df is None or len(reliability_df) == 0:
        logger.info(f"No data for reliability plot: {filename}")
        return

    if y_col not in reliability_df.columns:
        logger.warning(
            f"Cannot plot {filename}; missing column {y_col}"
        )
        return

    plot_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    plt.figure(figsize=(7, 5))

    plt.plot(
        reliability_df["mean_ood_percentile"],
        reliability_df[y_col],
        marker="o",
    )

    plt.xlabel(
        "OOD percentile (higher = less in-domain)"
    )

    plt.ylabel(ylabel)
    plt.title(title)

    plt.grid(True)
    plt.tight_layout()

    plt.savefig(
        plot_dir / filename,
        dpi=300,
    )

    plt.close()

def plot_reliability_curve_with_ood_boundaries(
    reliability_df,
    plot_dir,
    y_col,
    filename,
    ylabel,
    title,
    bin_boundaries=None,
):
    """
    Plot prediction reliability curve with vertical dotted OOD bin boundaries.

    X-axis:
        Mean OOD percentile.

    Y-axis:
        Mean error, e.g. mean_abs_error or mean_normalized_abs_error.

    Vertical dotted lines:
        OOD percentile boundaries used to define bins:
          in_domain / moderate_ood / high_ood / extreme_ood
    """

    if reliability_df is None or len(reliability_df) == 0:
        logger.info(f"No data for reliability plot: {filename}")
        return

    required_cols = {
        "mean_ood_percentile",
        y_col,
    }

    missing = required_cols - set(reliability_df.columns)

    if missing:
        logger.warning(
            f"Cannot plot {filename}; missing columns: {sorted(missing)}"
        )
        return

    if bin_boundaries is None:
        bin_boundaries = {
            "moderate_ood": 50.0,
            "high_ood": 80.0,
            "extreme_ood": 95.0,
        }

    plot_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    df = reliability_df.sort_values(
        "mean_ood_percentile"
    ).copy()

    fig, ax = plt.subplots(
        figsize=(7, 5)
    )

    ax.plot(
        df["mean_ood_percentile"],
        df[y_col],
        marker="o",
        color="#1D4ED8",
        linewidth=2,
        label=ylabel,
    )

    # -------------------------
    # Vertical dotted bin boundaries
    # -------------------------
    boundary_style = {
        "moderate_ood": {
            "color": "#7C3AED",
            "label": "moderate OOD boundary",
        },
        "high_ood": {
            "color": "#0F766E",
            "label": "high OOD boundary",
        },
        "extreme_ood": {
            "color": "#92400E",
            "label": "extreme OOD boundary",
        },
    }

    y_min, y_max = ax.get_ylim()
    y_text = y_min + 0.96 * (y_max - y_min)

    for boundary_name, boundary_value in bin_boundaries.items():
        style = boundary_style.get(
            boundary_name,
            {
                "color": "#6B7280",
                "label": boundary_name,
            },
        )

        ax.axvline(
            boundary_value,
            linestyle=":",
            linewidth=1.5,
            color=style["color"],
            alpha=0.8,
        )

        ax.text(
            boundary_value,
            y_text,
            f"{boundary_name}\n{boundary_value:.0f}%",
            rotation=90,
            va="top",
            ha="right",
            fontsize=8,
            color=style["color"],
        )

    ax.set_xlabel(
        "OOD percentile\n(higher = less in-domain)"
    )

    ax.set_ylabel(ylabel)
    ax.set_title(title)

    ax.set_xlim(0, 100)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    plt.savefig(
        plot_dir / filename,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()

def compute_task_ood_calibration_table(
    y_true,
    y_pred,
    task_names,
    ood_percentile,
    n_bins=10,
    min_n=5,
):
    """
    Per-task calibration of prediction error against OOD percentile.

    Unlike sampled reliability plots, this uses the full validation matrix.

    Outputs one row per task per OOD percentile bin.
    """

    rows = []

    ood_percentile = np.asarray(
        ood_percentile,
        dtype=float,
    )

    for task_idx, task in enumerate(task_names):
        yt = y_true[:, task_idx]
        yp = y_pred[:, task_idx]

        valid = (
            ~np.isnan(yt)
            & np.isfinite(ood_percentile)
        )

        if valid.sum() < min_n:
            continue

        task_df = pd.DataFrame({
            "ood_percentile": ood_percentile[valid],
            "y_true": yt[valid],
            "y_pred": yp[valid],
        })

        task_df["abs_error"] = np.abs(
            task_df["y_pred"] - task_df["y_true"]
        )

        task_df["ood_percentile_bin"] = pd.qcut(
            task_df["ood_percentile"],
            q=n_bins,
            duplicates="drop",
        )

        task_rows = []

        for b, sub in task_df.groupby(
            "ood_percentile_bin",
            observed=False,
        ):
            if len(sub) < min_n:
                continue

            metrics = compute_metrics(
                sub["y_true"].values,
                sub["y_pred"].values,
            )

            task_rows.append({
                "task": task,
                "ood_percentile_bin": str(b),
                "mean_ood_percentile": sub["ood_percentile"].mean(),
                "mean_abs_error": sub["abs_error"].mean(),
                "median_abs_error": sub["abs_error"].median(),
                "n_predictions": len(sub),
                **metrics,
            })

        if len(task_rows) == 0:
            continue

        task_calib = pd.DataFrame(task_rows)

        # Ratio to the lowest-OOD bin for this task.
        task_calib = task_calib.sort_values(
            "mean_ood_percentile"
        )

        baseline_mae = task_calib["mean_abs_error"].iloc[0]
        baseline_median_ae = task_calib["median_abs_error"].iloc[0]

        # Absolute degradation relative to the lowest-OOD bin.
        task_calib["mae_delta_vs_lowest_ood_bin"] = (
            task_calib["mean_abs_error"] - baseline_mae
        )

        task_calib["median_ae_delta_vs_lowest_ood_bin"] = (
            task_calib["median_abs_error"] - baseline_median_ae
        )

        # Relative degradation relative to the lowest-OOD bin.
        if np.isfinite(baseline_mae) and baseline_mae > 0:
            task_calib["mae_ratio_vs_lowest_ood_bin"] = (
                task_calib["mean_abs_error"] / baseline_mae
            )
        else:
            task_calib["mae_ratio_vs_lowest_ood_bin"] = np.nan

        if np.isfinite(baseline_median_ae) and baseline_median_ae > 0:
            task_calib["median_ae_ratio_vs_lowest_ood_bin"] = (
                task_calib["median_abs_error"] / baseline_median_ae
            )
        else:
            task_calib["median_ae_ratio_vs_lowest_ood_bin"] = np.nan

        rows.append(task_calib)

    if len(rows) == 0:
        return pd.DataFrame()

    return pd.concat(
        rows,
        ignore_index=True,
    )

def compute_task_reliability_curves(
    pred_long_df,
    n_bins=10,
    error_col="abs_error",
):
    rows = []

    required = {
        "task",
        "ood_percentile",
        error_col,
    }

    missing = required - set(pred_long_df.columns)

    if missing:
        logger.warning(
            f"Cannot compute task reliability curves; missing {sorted(missing)}"
        )
        return pd.DataFrame()

    for task, task_df in pred_long_df.groupby("task"):

        task_df = task_df[
            ["ood_percentile", error_col]
        ].dropna().copy()

        if len(task_df) < 5:
            continue

        task_df["confidence_bin"] = pd.qcut(
            task_df["ood_percentile"],
            q=n_bins,
            duplicates="drop",
        )

        for b, sub in task_df.groupby(
            "confidence_bin",
            observed=False,
        ):

            if len(sub) < 5:
                continue

            rows.append({
                "task": task,
                "confidence_bin": str(b),
                "mean_ood_percentile": sub["ood_percentile"].mean(),
                "mean_abs_error": sub[error_col].mean(),
                "n_predictions": len(sub),
                "error_col": error_col,
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

def plot_per_task_ood_degradation(
    task_ood_calibration_df,
    plot_dir,
    max_cols=3,
    ratio_thresholds=(1.25, 1.5, 2.0),
    delta_thresholds=None,
):
    """
    Plot OOD-driven degradation per task.

    Left y-axis:
        MAE ratio vs lowest-OOD bin.

    Right y-axis:
        MAE delta vs lowest-OOD bin.

    Dotted horizontal lines show interpretable thresholds.
    """

    if task_ood_calibration_df is None or len(task_ood_calibration_df) == 0:
        logger.info("No task OOD calibration data available for degradation plots")
        return

    required_cols = {
        "task",
        "mean_ood_percentile",
        "mae_ratio_vs_lowest_ood_bin",
        "mae_delta_vs_lowest_ood_bin",
    }

    missing = required_cols - set(task_ood_calibration_df.columns)

    if missing:
        logger.warning(
            "Cannot plot OOD degradation because columns are missing: "
            f"{sorted(missing)}"
        )
        return

    out_dir = plot_dir / "ood" / "ood_degradation"
    out_dir.mkdir(parents=True, exist_ok=True)

    tasks = list(task_ood_calibration_df["task"].dropna().unique())

    if len(tasks) == 0:
        return

    # -------------------------
    # Individual task plots
    # -------------------------
    for task in tasks:
        task_df = (
            task_ood_calibration_df[
                task_ood_calibration_df["task"] == task
            ]
            .sort_values("mean_ood_percentile")
            .copy()
        )

        if len(task_df) == 0:
            continue

        fig, ax1 = plt.subplots(figsize=(7, 5))

        ax1.plot(
            task_df["mean_ood_percentile"],
            task_df["mae_ratio_vs_lowest_ood_bin"],
            marker="o",
            color="#1D4ED8",
            label="MAE ratio",
        )

        ax1.set_xlabel("Mean OOD percentile")
        ax1.set_ylabel("MAE ratio vs lowest-OOD bin", color="#1D4ED8")
        ax1.tick_params(axis="y", labelcolor="#1D4ED8")
        ax1.grid(True, alpha=0.3)

        for thr in ratio_thresholds:
            ax1.axhline(
                thr,
                linestyle=":",
                linewidth=1.2,
                color="#1D4ED8",
                alpha=0.6,
            )

            ax1.text(
                task_df["mean_ood_percentile"].min(),
                thr,
                f"{thr:.2f}x",
                color="#1D4ED8",
                fontsize=8,
                va="bottom",
            )

        ax2 = ax1.twinx()

        ax2.plot(
            task_df["mean_ood_percentile"],
            task_df["mae_delta_vs_lowest_ood_bin"],
            marker="s",
            color="#0F766E",
            label="MAE delta",
        )

        ax2.set_ylabel("MAE delta vs lowest-OOD bin", color="#0F766E")
        ax2.tick_params(axis="y", labelcolor="#0F766E")

        if delta_thresholds is not None:
            for thr in delta_thresholds:
                ax2.axhline(
                    thr,
                    linestyle=":",
                    linewidth=1.2,
                    color="#0F766E",
                    alpha=0.6,
                )

                ax2.text(
                    task_df["mean_ood_percentile"].max(),
                    thr,
                    f"+{thr:.2f}",
                    color="#0F766E",
                    fontsize=8,
                    va="bottom",
                    ha="right",
                )

        lines_1, labels_1 = ax1.get_legend_handles_labels()
        lines_2, labels_2 = ax2.get_legend_handles_labels()

        ax1.legend(
            lines_1 + lines_2,
            labels_1 + labels_2,
            loc="upper left",
            frameon=False,
        )

        plt.title(task)
        plt.tight_layout()

        safe_task = (
            str(task)
            .replace("/", "_")
            .replace("\\", "_")
            .replace(" ", "_")
            .replace(":", "_")
        )

        plt.savefig(
            out_dir / f"{safe_task}_ood_degradation.png",
            dpi=300,
            bbox_inches="tight",
        )

        plt.close()

    # -------------------------
    # Combined grid plot
    # -------------------------
    n_tasks = len(tasks)
    n_rows = math.ceil(n_tasks / max_cols)

    fig, axes = plt.subplots(
        n_rows,
        max_cols,
        figsize=(5 * max_cols, 4 * n_rows),
    )

    axes = np.array(axes).flatten()

    for i, task in enumerate(tasks):
        ax = axes[i]

        task_df = (
            task_ood_calibration_df[
                task_ood_calibration_df["task"] == task
            ]
            .sort_values("mean_ood_percentile")
            .copy()
        )

        ax.plot(
            task_df["mean_ood_percentile"],
            task_df["mae_ratio_vs_lowest_ood_bin"],
            marker="o",
            color="#1D4ED8",
            label="MAE ratio",
        )

        for thr in ratio_thresholds:
            ax.axhline(
                thr,
                linestyle=":",
                linewidth=1.0,
                color="#1D4ED8",
                alpha=0.5,
            )

        ax.set_title(task)
        ax.set_xlabel("OOD percentile")
        ax.set_ylabel("MAE ratio")
        ax.grid(True, alpha=0.3)

    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()

    plt.savefig(
        out_dir / "per_task_ood_mae_ratio_grid.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()

def summarise_task_ood_degradation(
    task_ood_calibration_df,
    ratio_thresholds=(1.25, 1.5, 2.0),
    delta_thresholds=None,
):
    """
    Summarise OOD degradation per task.

    Reports:
      - maximum MAE ratio
      - maximum MAE delta
      - first OOD percentile crossing each ratio threshold
      - first OOD percentile crossing each delta threshold, if provided
    """

    if task_ood_calibration_df is None or len(task_ood_calibration_df) == 0:
        return pd.DataFrame()

    rows = []

    for task, task_df in task_ood_calibration_df.groupby("task"):
        task_df = task_df.sort_values("mean_ood_percentile").copy()

        row = {
            "task": task,
            "n_bins": len(task_df),
            "max_mae_ratio_vs_lowest_ood_bin": task_df[
                "mae_ratio_vs_lowest_ood_bin"
            ].max(),
            "max_mae_delta_vs_lowest_ood_bin": task_df[
                "mae_delta_vs_lowest_ood_bin"
            ].max(),
            "mean_mae_ratio_vs_lowest_ood_bin": task_df[
                "mae_ratio_vs_lowest_ood_bin"
            ].mean(),
            "mean_mae_delta_vs_lowest_ood_bin": task_df[
                "mae_delta_vs_lowest_ood_bin"
            ].mean(),
        }

        for thr in ratio_thresholds:
            crossing = task_df[
                task_df["mae_ratio_vs_lowest_ood_bin"] >= thr
            ]

            col = f"first_ood_percentile_mae_ratio_ge_{thr}"

            if len(crossing) > 0:
                row[col] = crossing["mean_ood_percentile"].iloc[0]
            else:
                row[col] = np.nan

        if delta_thresholds is not None:
            for thr in delta_thresholds:
                crossing = task_df[
                    task_df["mae_delta_vs_lowest_ood_bin"] >= thr
                ]

                col = f"first_ood_percentile_mae_delta_ge_{thr}"

                if len(crossing) > 0:
                    row[col] = crossing["mean_ood_percentile"].iloc[0]
                else:
                    row[col] = np.nan

        rows.append(row)

    return pd.DataFrame(rows)

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

    logger.info("Computing task metrics")

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

    logger.info("Plotting task correlation heatmap")

    max_corr_rows = int(config.get("max_corr_rows", 50000))

    if y_true.shape[0] > max_corr_rows:
        rng = np.random.default_rng(42)
        idx = rng.choice(
            y_true.shape[0],
            size=max_corr_rows,
            replace=False,
        )
        corr_data = y_true[idx, :]
    else:
        corr_data = y_true

    df_tasks = pd.DataFrame(
        corr_data,
        columns=task_names,
    )

    plt.figure(figsize=(6, 5))
    sns.heatmap(
        df_tasks.corr(),
        annot=True,
        cmap="coolwarm_r",
        center=0,
    )
    plt.title("Task correlation matrix")
    plt.tight_layout()
    plt.savefig(plot_dir / "task_correlation.png", dpi=300)
    plt.close()

    # =========================
    # 3. TRUE VS PREDICTED
    # =========================

    logger.info("Plotting true vs predicted panels")

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

    logger.info("Plotting residual panels")

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

    logger.info("Plotting residual histograms")

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

    logger.info("Starting OOD evaluation block")

    ood_cfg = config.get("ood", {})

    if ood_cfg.get("enabled", False):
        logger.info("Computing OOD/applicability-domain metrics")

        train_smiles, val_smiles = load_split_smiles_from_context_or_disk(
            context=context,
            config=config,
        )

        if len(val_smiles) != y_true.shape[0]:
            raise ValueError(
                "Number of n SMILES does not match y_true rows: "
                f"len(val_smiles)={len(val_smiles)}, "
                f"y_true.shape[0]={y_true.shape[0]}. "
                "Check that validation predictions were generated with "
                "shuffle=False and that the SMILES metadata corresponds "
                "to this prediction run."
            )

        # --------------------------------------------------
        # Molecule-level OOD features only.
        # One row per validation molecule.
        # --------------------------------------------------
        ood_df = compute_ood_features(
            query_smiles=val_smiles,
            train_smiles=train_smiles,
            config=ood_cfg,
        )

        ood_df.insert(
            0,
            "smiles",
            val_smiles,
        )

        ood_df.to_csv(
            plot_dir / "ood_features_by_molecule.csv",
            index=False,
        )

        logger.info(
            f"Wrote molecule-level OOD features: "
            f"{plot_dir / 'ood_features_by_molecule.csv'}"
        )

        # --------------------------------------------------
        # Metrics by OOD bin, no long dataframe.
        # --------------------------------------------------
        ood_bin_values = ood_df["ood_bin"].values

        per_task_ood = metrics_by_molecule_group(
            y_true=y_true,
            y_pred=y_pred,
            task_names=task_names,
            group_values=ood_bin_values,
            group_name="ood_bin",
            min_n=5,
        )

        overall_ood = overall_metrics_by_molecule_group(
            y_true=y_true,
            y_pred=y_pred,
            group_values=ood_bin_values,
            group_name="ood_bin",
            min_n=20,
        )

        metrics_by_ood_bin = pd.concat(
            [overall_ood, per_task_ood],
            ignore_index=True,
        )

        metrics_by_ood_bin.to_csv(
            plot_dir / "metrics_by_ood_bin.csv",
            index=False,
        )

        # --------------------------------------------------
        # OOD component summary.
        # Helps interpret what drives the composite score.
        # --------------------------------------------------
        component_cols = [
            c for c in ood_df.columns
            if c.startswith("ood_component_")
        ]

        if len(component_cols) > 0:
            component_summary = (
                ood_df[component_cols + ["ood_score", "ood_percentile"]]
                .describe()
                .T
                .reset_index()
                .rename(columns={"index": "feature"})
            )

            component_summary.to_csv(
                plot_dir / "ood_component_summary.csv",
                index=False,
            )

        # --------------------------------------------------
        # Known vs novel scaffold.
        # --------------------------------------------------
        scaffold_metrics_per_task = metrics_by_molecule_group(
            y_true=y_true,
            y_pred=y_pred,
            task_names=task_names,
            group_values=ood_df["known_scaffold"].values,
            group_name="known_scaffold",
            min_n=5,
        )

        scaffold_metrics_overall = overall_metrics_by_molecule_group(
            y_true=y_true,
            y_pred=y_pred,
            group_values=ood_df["known_scaffold"].values,
            group_name="known_scaffold",
            min_n=20,
        )

        scaffold_metrics = pd.concat(
            [scaffold_metrics_overall, scaffold_metrics_per_task],
            ignore_index=True,
        )

        scaffold_metrics.to_csv(
            plot_dir / "metrics_by_scaffold.csv",
            index=False,
        )

        # --------------------------------------------------
        # Similarity bins.
        # --------------------------------------------------
        similarity_metrics = metrics_by_continuous_molecule_feature(
            y_true=y_true,
            y_pred=y_pred,
            task_names=task_names,
            values=ood_df["nearest_train_tanimoto"].values,
            feature_name="similarity",
            bins=np.linspace(0.0, 1.0, 6),
            labels=[
                "0.0-0.2",
                "0.2-0.4",
                "0.4-0.6",
                "0.6-0.8",
                "0.8-1.0",
            ],
            min_n=5,
        )

        similarity_metrics.to_csv(
            plot_dir / "metrics_by_similarity.csv",
            index=False,
        )

        # --------------------------------------------------
        # Physchem-distance quantile bins.
        # --------------------------------------------------
        distance_metrics = metrics_by_quantile_molecule_feature(
            y_true=y_true,
            y_pred=y_pred,
            task_names=task_names,
            values=ood_df["physchem_robust_distance"].values,
            feature_name="physchem_distance",
            q=5,
            min_n=5,
        )

        distance_metrics.to_csv(
            plot_dir / "metrics_by_physchem_distance.csv",
            index=False,
        )

        # --------------------------------------------------
        # Local density bins using neighbour counts.
        # --------------------------------------------------
        density_col = ood_cfg.get(
            "density_metric_col",
            "n_train_neighbors_ge_0.6",
        )

        if density_col in ood_df.columns:
            density_metrics = metrics_by_quantile_molecule_feature(
                y_true=y_true,
                y_pred=y_pred,
                task_names=task_names,
                values=ood_df[density_col].values,
                feature_name=density_col,
                q=5,
                min_n=5,
            )

            density_metrics.to_csv(
                plot_dir / "metrics_by_local_density.csv",
                index=False,
            )

        # --------------------------------------------------
        # Sampled long dataframe for plotting only.
        # --------------------------------------------------
        sample_n = int(
            ood_cfg.get(
                "plot_sample_rows",
                100000,
            )
        )

        task_error_scale_method = ood_cfg.get(
            "task_error_scale_method",
            "mae",
        )

        task_error_scales = compute_task_error_scales(
            y_true=y_true,
            y_pred=y_pred,
            task_names=task_names,
            method=task_error_scale_method,
        )

        pd.DataFrame([
            {
                "task": task,
                "error_scale_method": task_error_scale_method,
                "error_scale": scale,
            }
            for task, scale in task_error_scales.items()
        ]).to_csv(
            plot_dir / "task_error_scales.csv",
            index=False,
        )

        pred_sample_df = build_prediction_sample_table(
            y_true=y_true,
            y_pred=y_pred,
            task_names=task_names,
            val_smiles=val_smiles,
            ood_df=ood_df,
            task_error_scales=task_error_scales,
            max_rows=sample_n,
            random_seed=42,
        )

        pred_sample_df.to_csv(
            plot_dir / "predictions_with_ood_scores_sampled_long.csv",
            index=False,
        )

        task_ood_calibration = compute_task_ood_calibration_table(
            y_true=y_true,
            y_pred=y_pred,
            task_names=task_names,
            ood_percentile=ood_df["ood_percentile"].values,
            n_bins=ood_cfg.get("task_calibration_bins", 10),
            min_n=ood_cfg.get("task_calibration_min_n", 5),
        )

        task_ood_calibration.to_csv(
            plot_dir / "task_ood_calibration.csv",
            index=False,
        )

        if len(task_ood_calibration) > 0:
            plot_per_task_reliability_curves(
                task_reliability_df=task_ood_calibration,
                plot_dir=plot_dir,
            )

            ratio_thresholds = tuple(
                ood_cfg.get(
                    "mae_ratio_thresholds",
                    [1.25, 1.5, 2.0],
                )
            )

            delta_thresholds_cfg = ood_cfg.get(
                "mae_delta_thresholds",
                None,
            )

            if delta_thresholds_cfg is not None:
                delta_thresholds = tuple(delta_thresholds_cfg)
            else:
                delta_thresholds = None

            plot_per_task_ood_degradation(
                task_ood_calibration_df=task_ood_calibration,
                plot_dir=plot_dir,
                max_cols=3,
                ratio_thresholds=ratio_thresholds,
                delta_thresholds=delta_thresholds,
            )

            task_ood_degradation_summary = summarise_task_ood_degradation(
                task_ood_calibration_df=task_ood_calibration,
                ratio_thresholds=ratio_thresholds,
                delta_thresholds=delta_thresholds,
            )

            task_ood_degradation_summary.to_csv(
                plot_dir / "task_ood_degradation_summary.csv",
                index=False,
            )

        if len(pred_sample_df) > 0:
            density_col = ood_cfg.get(
                "density_metric_col",
                "n_train_neighbors_ge_0.6",
            )

            corr_df = compute_error_correlations(
                pred_sample_df,
                candidate_features=[
                    "nearest_train_tanimoto",
                    "mean_top5_train_tanimoto",
                    "physchem_robust_distance",
                    density_col,
                    "ood_score",
                    "ood_percentile",
                ],
            )


            corr_df.to_csv(
                plot_dir / "ood_feature_correlations_sampled.csv",
                index=False,
            )

            reliability_df = compute_reliability_curve(
                pred_sample_df,
                confidence_col="ood_percentile",
                error_col="abs_error",
                n_bins=10,
            )

            reliability_df.to_csv(
                plot_dir / "prediction_reliability_curve_sampled.csv",
                index=False,
            )

            normalized_reliability_df = compute_reliability_curve(
                pred_sample_df,
                confidence_col="ood_percentile",
                error_col="normalized_abs_error",
                n_bins=10,
            )

            normalized_reliability_df.to_csv(
                plot_dir / "prediction_reliability_curve_normalized_sampled.csv",
                index=False,
            )

            task_reliability = compute_task_reliability_curves(
                pred_sample_df,
                n_bins=10,
                error_col="abs_error",
            )

            task_reliability.to_csv(
                plot_dir / "task_prediction_reliability_curves_sampled.csv",
                index=False,
            )

            task_reliability_norm = compute_task_reliability_curves(
                pred_sample_df,
                n_bins=10,
                error_col="normalized_abs_error",
            )

            task_reliability_norm.to_csv(
                plot_dir / "task_prediction_reliability_curves_normalized_sampled.csv",
                index=False,
            )

            plot_ood_diagnostics(
                pred_long_df=pred_sample_df,
                plot_dir=plot_dir,
            )

            plot_per_task_ood_panels(
                pred_long_df=pred_sample_df,
                plot_dir=plot_dir,
                task_names=task_names,
            )

            bin_cfg = ood_cfg.get("binning", {})

            ood_bin_boundaries = {
                "moderate_ood": bin_cfg.get("moderate_cut", 50.0),
                "high_ood": bin_cfg.get("high_cut", 80.0),
                "extreme_ood": bin_cfg.get("extreme_cut", 95.0),
            }

            plot_reliability_curve_with_ood_boundaries(
                reliability_df,
                plot_dir / "ood",
                y_col="mean_abs_error",
                filename="prediction_reliability_curve.png",
                ylabel="Mean absolute error",
                title="Prediction reliability curve",
                bin_boundaries=ood_bin_boundaries,
            )

            plot_reliability_curve_with_ood_boundaries(
                normalized_reliability_df,
                plot_dir / "ood",
                y_col="mean_normalized_abs_error",
                filename="prediction_reliability_curve_normalized.png",
                ylabel="Mean normalized absolute error",
                title="Prediction reliability curve, task-normalized",
                bin_boundaries=ood_bin_boundaries,
            )

            plot_reliability_curve_generic(
                normalized_reliability_df,
                plot_dir / "ood",
                y_col="mean_normalized_abs_error",
                filename="prediction_reliability_curve_normalized_plain.png",
                ylabel="Mean normalized absolute error",
                title="Prediction reliability curve, task-normalized",
            )

            plot_error_correlations(
                corr_df,
                plot_dir / "ood",
            )

            plot_per_task_reliability_curves(
                task_reliability_df=task_reliability,
                plot_dir=plot_dir,
            )

        logger.info(f"OOD evaluation written to {plot_dir / 'ood'}")

    logger.info(f"Saved evaluation plots to {plot_dir}")

    return {}