import yaml
import re
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import root_mean_squared_error, mean_absolute_error, r2_score
import matplotlib.pyplot as plt
from dataclasses import dataclass

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@dataclass
class BootstrapMetric:
    mean: float
    low: float
    high: float

@dataclass
class EvaluationResult:
    metrics: dict
    predictions: pd.DataFrame

def safe_filename(s: str) -> str:
    """
    Convert arbitrary strings to filesystem-safe filenames.
    """
    s = s.replace("/", "_per_")
    s = re.sub(r"[^\w\d\-_. ]+", "", s)
    s = s.replace(" ", "_")
    return s

def _bootstrap_metric(y_true, y_pred, fn, n_boot=500, seed=42):
    rng = np.random.default_rng(seed)
    stats = []

    idx = np.arange(len(y_true))
    for _ in range(n_boot):
        sample = rng.choice(idx, size=len(idx), replace=True)
        stats.append(fn(y_true[sample], y_pred[sample]))

    stats = np.array(stats)
    return BootstrapMetric(
        mean=float(stats.mean()),
        low=float(np.percentile(stats, 2.5)),
        high=float(np.percentile(stats, 97.5))
    )


def evaluate_with_bootstrap(
    models,
    pca,
    embed_fn,
    df_eval,
    properties,
    n_boot=500,
    min_n=20
):
    # --- Embed once on the correct device ---
    X = embed_fn(df_eval["sequence"].tolist())

    if pca is not None:
        X_pca = pca.transform(X)
    else:
        # Already PCA-transformed
        X_pca = X

    # Precompute predictions for all sequences
    preds = {
        prop: model.predict(X_pca)
        for prop, model in models.items()
    }

    df_pred = df_eval.copy()
    for prop in properties:
        df_pred[f"{prop}_pred"] = preds[prop]

    metrics = {}

    for prop in properties:
        y_true = df_pred[prop].values
        y_pred = df_pred[f"{prop}_pred"].values
        mask = np.isfinite(y_true)

        if mask.sum() < 2:
            logger.warning(f"Skipping {prop}: not enough data to fit even minimal model")
            continue

        yt = y_true[mask]
        yp = y_pred[mask]

        # --- Bootstrapping using numpy (CPU only) ---
        rmse = _bootstrap_metric(
            yt, yp,
            lambda a, b: np.sqrt(np.mean((a - b) ** 2)),
            n_boot=n_boot
        )

        mae = _bootstrap_metric(
            yt, yp,
            lambda a, b: np.mean(np.abs(a - b)),
            n_boot=n_boot
        )

        bias = _bootstrap_metric(
            yt, yp,
            lambda a, b: np.mean(b - a),
            n_boot=n_boot
        )

        pearson = _bootstrap_metric(
            yt, yp,
            lambda a, b: pearsonr(a, b)[0],
            n_boot=n_boot
        )

        spearman = _bootstrap_metric(
            yt, yp,
            lambda a, b: spearmanr(a, b)[0],
            n_boot=n_boot
        )

        # Calibration (point estimates only)
        slope, intercept = np.polyfit(yp, yt, 1)
        r2 = np.corrcoef(yp, yt)[0, 1] ** 2  # simple R²

        metrics[prop] = {
            "rmse": rmse,
            "mae": mae,
            "bias": bias,
            "pearson": pearson,
            "spearman": spearman,
            "calibration": {
                "slope": float(slope),
                "intercept": float(intercept),
                "r2": float(r2),
            },
            "n": int(mask.sum()),
        }

    return EvaluationResult(metrics=metrics, predictions=df_pred)



def save_evaluation(result, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    plots_dir = out_dir / "plots"
    plots_dir.mkdir(exist_ok=True)

    # --- Metrics ---
    serialisable = {
        prop: {
            k: vars(v) if hasattr(v, "__dict__") else v
            for k, v in metrics.items()
        }
        for prop, metrics in result.metrics.items()
    }

    with open(out_dir / "metrics.yaml", "w") as f:
        yaml.safe_dump(serialisable, f)

    # --- Predictions ---
    df_full = result.predictions.copy()

    # --- Fill experimental NaNs with predictions ---
    for prop in result.metrics.keys():
        exp_col = prop
        pred_col = f"{prop}_pred"
        if exp_col in df_full.columns and pred_col in df_full.columns:
            mask = ~np.isfinite(df_full[exp_col])
            df_full.loc[mask, exp_col] = df_full.loc[mask, pred_col]

    # Save full prediction surface (ALL rows × ALL predicted properties)
    df_full.to_csv(out_dir / "all_predictions_all_properties.csv", index=False)

    # Save masked predictions (only rows with at least one experimental value)
    mask_any = np.zeros(len(df_full), dtype=bool)
    for col in df_full.columns:
        if not col.endswith("_pred") and col not in ["sequence", "ancestor_id"]:
            mask_any |= np.isfinite(df_full[col].values)
    df_masked = df_full.loc[mask_any].copy()
    df_masked.to_csv(out_dir / "predictions.csv", index=False)

    # --- Plots ---
    for prop, metrics in result.metrics.items():
        if prop not in df_full:
            continue

        y_true = df_full[prop].values
        y_pred = df_full[f"{prop}_pred"].values
        mask = np.isfinite(y_true)

        if mask.sum() < 5:
            continue

        yt = y_true[mask]
        yp = y_pred[mask]

        fname = safe_filename(prop)

        # Parity plot
        plot_parity(
            yt,
            yp,
            prop,
            metrics["rmse"],
            plots_dir / f"{fname}_parity.png"
        )

        # Error plot
        plot_error_vs_value(
            yt,
            yp,
            prop,
            plots_dir / f"{fname}_error.png"
        )


def plot_parity(
    y_true,
    y_pred,
    prop,
    rmse_metric,
    out_path
):
    rmse = rmse_metric.mean

    plt.figure(figsize=(5, 5))
    plt.scatter(y_true, y_pred, alpha=0.6)

    lims = [
        min(y_true.min(), y_pred.min()),
        max(y_true.max(), y_pred.max())
    ]
    plt.plot(lims, lims, "k--", alpha=0.8, label="Ideal")

    plt.fill_between(
        lims,
        [l + rmse for l in lims],
        [l - rmse for l in lims],
        color="gray",
        alpha=0.15,
        label="±1 RMSE"
    )

    plt.xlabel(f"Experimental {prop}")
    plt.ylabel(f"Predicted {prop}")
    plt.title(f"{prop}: predicted vs experimental")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


def plot_error_vs_value(
    y_true,
    y_pred,
    prop,
    out_path
):
    err = y_pred - y_true

    plt.figure(figsize=(5, 4))
    plt.scatter(y_true, err, alpha=0.6)
    plt.axhline(0, color="k", linestyle="--", alpha=0.8)

    plt.xlabel(f"Experimental {prop}")
    plt.ylabel("Prediction error")
    plt.title(f"{prop}: error vs value")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()
