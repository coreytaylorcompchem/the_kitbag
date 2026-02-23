import yaml
import re
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import root_mean_squared_error, mean_absolute_error, r2_score
import matplotlib.pyplot as plt
from dataclasses import dataclass

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def debug_numeric_array(name, arr, max_examples=5):
    logger.info(f"[DEBUG] {name}: dtype={arr.dtype}, shape={arr.shape}")

    # Check if object dtype
    if arr.dtype == object:
        logger.error(f"[DEBUG] {name} is OBJECT dtype")

        # Find first few non-numeric entries
        bad = []
        for i, v in enumerate(arr):
            if not isinstance(v, (int, float, np.integer, np.floating)) and not pd.isna(v):
                bad.append((i, type(v), v))
            if len(bad) >= max_examples:
                break

        logger.error(f"[DEBUG] {name} sample non-numeric entries: {bad}")

    # Test np.isfinite safely
    try:
        mask = np.isfinite(arr)
        logger.info(f"[DEBUG] {name} np.isfinite OK, finite_count={mask.sum()}")
    except Exception as e:
        logger.error(f"[DEBUG] np.isfinite FAILED on {name}: {e}")

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
    preds = {}

    df_pred = df_eval.copy()

    for prop in properties:

        model_list = models[prop]

        all_preds = np.column_stack([
            m.predict(X_pca) for m in model_list
        ])

        mean_pred = all_preds.mean(axis=1)
        std_pred  = all_preds.std(axis=1)

        df_pred[f"{prop}_pred"]     = mean_pred
        df_pred[f"{prop}_pred_std"] = std_pred
        df_pred[f"{prop}_pred_lci"] = mean_pred - 1.96 * std_pred
        df_pred[f"{prop}_pred_uci"] = mean_pred + 1.96 * std_pred

        preds[prop] = mean_pred

    metrics = {}

    for prop in properties:
        y_true = df_pred[prop].values
        y_pred = df_pred[f"{prop}_pred"].values

        logger.info(f"[DEBUG] Eval property {prop}")
        debug_numeric_array(f"eval_{prop}_y_true", y_true)

        try:
            mask = np.isfinite(y_true)
        except Exception as e:
            logger.error(f"[FATAL DEBUG] np.isfinite crashed in evaluation for {prop}")
            debug_numeric_array(f"eval_{prop}_CRASH", y_true)
            raise

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



def save_evaluation(result, out_dir, score_config=None, train_stats=None):
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

    # ----------------------------------------------------------
    # Composite scoring (optional)
    # ----------------------------------------------------------
    if score_config is not None and train_stats is not None:

        method = score_config.get("method", "simple")
        mutation_penalty = score_config.get("mutation_penalty", 0.0)
        prop_weights = score_config.get("property_weights", {})
        opt_direction = score_config.get("optimisation_direction", {})
        lower_is_better = set(score_config.get("lower_is_better", []))

        scores = np.zeros(len(df_full))

        if method == "zscore_weighted":
            for prop in result.metrics.keys():
                pred_col = f"{prop}_pred"

                if pred_col not in df_full.columns:
                    continue

                if prop not in train_stats.columns:
                    continue

                mu = train_stats.loc["mean", prop]
                sigma = train_stats.loc["std", prop] + 1e-8

                z = (df_full[pred_col] - mu) / sigma

                # flip if lower is better
                direction = opt_direction.get(prop, 1)
                if prop in lower_is_better:
                    direction *= -1

                weight = prop_weights.get(prop, 1.0)
                scores += weight * direction * z

            # apply mutation penalty if available
            scores -= mutation_penalty * df_full.get("mut_count", 0)

        elif method == "simple":
            # fallback: single-property simple score
            scores = df_full.get("Tm1 (°C)", 0) - mutation_penalty * df_full.get("mut_count", 0)

        df_full["composite_score"] = scores  # add column before saving
    
    if (
        score_config is not None
        and train_stats is not None
        and "composite_score" in df_full.columns
    ):

        method = score_config.get("method", "simple")
        prop_weights = score_config.get("property_weights", {})
        opt_direction = score_config.get("optimisation_direction", {})
        lower_is_better = set(score_config.get("lower_is_better", []))

        score_var = np.zeros(len(df_full))

        # ----------------------------------------------------------
        # Z-SCORED MULTI-PROPERTY
        # ----------------------------------------------------------
        if method == "zscore_weighted":

            for prop in result.metrics.keys():

                std_col = f"{prop}_pred_std"
                pred_col = f"{prop}_pred"

                if std_col not in df_full or prop not in train_stats.columns:
                    continue

                mu = train_stats.loc["mean", prop]
                sigma = train_stats.loc["std", prop] + 1e-8

                weight = prop_weights.get(prop, 1.0)
                direction = opt_direction.get(prop, 1)

                if prop in lower_is_better:
                    direction *= -1

                # Variance propagation:
                # composite = Σ w * direction * (pred - mu)/sigma
                # Var = Σ (w*direction/sigma)^2 * Var(pred)

                coeff = weight * direction / sigma
                score_var += (coeff ** 2) * (df_full[std_col] ** 2)

        # ----------------------------------------------------------
        # SIMPLE SCORING
        # ----------------------------------------------------------
        elif method == "simple":

            # Example: Tm1 - penalty * mut_count
            # So variance just comes from the predicted property

            main_prop = "Tm1 (°C)"

            std_col = f"{main_prop}_pred_std"

            if std_col in df_full.columns:
                score_var = df_full[std_col] ** 2
            else:
                score_var = np.zeros(len(df_full))

        # ----------------------------------------------------------
        # PARETO RANK
        # ----------------------------------------------------------
        elif method == "pareto_rank":

            # Rank-based score → not analytically differentiable
            # Uncertainty propagation not meaningful
            # Set to NaN or zero

            score_var = np.full(len(df_full), np.nan)

        # ----------------------------------------------------------
        # Finalise CI
        # ----------------------------------------------------------

        score_std = np.sqrt(score_var)

        df_full["composite_score_std"] = score_std
        df_full["composite_score_lci"] = df_full["composite_score"] - 1.96 * score_std
        df_full["composite_score_uci"] = df_full["composite_score"] + 1.96 * score_std

    # --- Prepare the all-properties CSV ---
    df_all_props = df_full.copy()
    # Fill experimental NaNs with predictions only for the all-properties CSV
    for prop in result.metrics.keys():
        exp_col = prop
        pred_col = f"{prop}_pred"
        if exp_col in df_all_props.columns and pred_col in df_all_props.columns:
            if np.issubdtype(df_all_props[exp_col].dtype, np.number):
                mask = ~np.isfinite(df_all_props[exp_col].values)
                df_all_props.loc[mask, exp_col] = df_all_props.loc[mask, pred_col]
            df_all_props.loc[mask, exp_col] = df_all_props.loc[mask, pred_col]

    # Save full prediction surface (ALL rows × ALL predicted properties)
    df_all_props.to_csv(out_dir / "all_predictions_all_properties.csv", index=False)

    # --- Prepare the masked predictions CSV ---
    mask_any = np.zeros(len(df_full), dtype=bool)

    for col in df_full.columns:

        # Skip non-numeric columns entirely
        if not np.issubdtype(df_full[col].dtype, np.number):
            continue

        if not col.endswith("_pred"):
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
            plots_dir / f"{fname}_parity.png",
            df_full=df_full.loc[mask]
        )

        # Error plot
        plot_error_vs_value(
            yt,
            yp,
            prop,
            plots_dir / f"{fname}_error.png",
            df_subset=df_full.loc[mask]
        )



def plot_parity(
    y_true,
    y_pred,
    prop,
    rmse_metric,
    out_path,
    df_subset=None   # already masked rows
):
    rmse = rmse_metric.mean

    plt.figure(figsize=(5, 5))

    # --- Source-coloured plotting ---
    if df_subset is not None and "source" in df_subset.columns:

        sources = df_subset["source"].values
        unique_sources = sorted(pd.unique(sources))

        palette = sns.color_palette("tab10", len(unique_sources))
        colour_map = dict(zip(unique_sources, palette))

        for src in unique_sources:
            idx = sources == src
            plt.scatter(
                y_true[idx],
                y_pred[idx],
                label=str(src),
                alpha=0.7
            )

        plt.legend(title="Source")

    else:
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
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


def plot_error_vs_value(
    y_true,
    y_pred,
    prop,
    out_path,
    df_subset=None
):
    err = y_pred - y_true

    plt.figure(figsize=(5, 4))

    if df_subset is not None and "source" in df_subset.columns:

        sources = df_subset["source"].values
        unique_sources = sorted(pd.unique(sources))

        palette = sns.color_palette("tab10", len(unique_sources))
        colour_map = dict(zip(unique_sources, palette))

        for src in unique_sources:
            idx = sources == src
            plt.scatter(
                y_true[idx],
                err[idx],
                label=str(src),
                alpha=0.7
            )

        plt.legend(title="Source")

    else:
        plt.scatter(y_true, err, alpha=0.6)

    plt.axhline(0, color="k", linestyle="--", alpha=0.8)

    plt.xlabel(f"Experimental {prop}")
    plt.ylabel("Prediction error")
    plt.title(f"{prop}: error vs value")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()
