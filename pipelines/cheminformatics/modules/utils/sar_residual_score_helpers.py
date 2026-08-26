import os

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

from joblib import Parallel, delayed, parallel_config

from rdkit.Chem import Draw
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw, Descriptors, Crippen, rdMolDescriptors, DataStructs

from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor, HistGradientBoostingRegressor
from sklearn.model_selection import KFold
from sklearn.base import clone
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error

from scipy.stats import norm

from modules.utils.free_wilson_helpers import get_rgroup_columns

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def plot_virtual_prediction_components(
    virtual_df,
    output_path,
    top_n=100,
    sort_col="rank_predicted_potency",
    dpi=300,
):
    """
    Plot FW baseline, ML residual correction, and combined prediction
    for top virtual candidates.
    """

    if virtual_df is None or virtual_df.empty:
        return None

    required = {
        "fw_predicted_activity",
        "ml_residual_pred",
        "combined_predicted_activity",
        sort_col,
    }

    missing = required - set(virtual_df.columns)
    if missing:
        logger.warning(
            f"[FW residual ML] Component plot skipped, missing columns: {missing}"
        )
        return None

    df = virtual_df.copy()

    for col in ["fw_predicted_activity", "ml_residual_pred", "combined_predicted_activity", sort_col]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["fw_predicted_activity", "ml_residual_pred", "combined_predicted_activity", sort_col])

    if df.empty:
        return None

    sort_ascending = str(sort_col).startswith("rank_")

    df = (
        df.sort_values(
            sort_col,
            ascending=sort_ascending,
        )
        .head(top_n)
        .reset_index(drop=True)
    )

    df["rank"] = np.arange(1, len(df) + 1)

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    axes[0].plot(df["rank"], df["fw_predicted_activity"], label="FW baseline", linewidth=2)
    axes[0].plot(df["rank"], df["combined_predicted_activity"], label="FW + ML residual", linewidth=2)
    axes[0].set_ylabel("Predicted activity")
    axes[0].set_title(f"Top {len(df)} virtual candidates: prediction components")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    colors = ["#2ca25f" if v >= 0 else "#de2d26" for v in df["ml_residual_pred"]]

    axes[1].bar(df["rank"], df["ml_residual_pred"], color=colors)
    axes[1].axhline(0, color="black", linewidth=1)
    axes[1].set_xlabel(f"Candidates ordered by {sort_col}")
    axes[1].set_ylabel("ML residual correction")
    axes[1].set_title("ML residual contribution")
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path

def plot_observed_fw_vs_ml_validation(
    observed_df,
    output_path,
    activity_col,
    fw_pred_col="fw_predicted_activity",
    combined_pred_col="combined_predicted_activity_cv",
    dpi=300,
):
    """
    Compare FW-only and FW+ML predictions on observed compounds.
    Uses CV combined predictions if available.
    """

    if observed_df is None or observed_df.empty:
        return None

    required = {activity_col, fw_pred_col, combined_pred_col}
    missing = required - set(observed_df.columns)

    if missing:
        logger.warning(
            f"[FW residual ML] Validation plot skipped, missing columns: {missing}"
        )
        return None

    df = observed_df.copy()

    for col in [activity_col, fw_pred_col, combined_pred_col]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=[activity_col, fw_pred_col, combined_pred_col])

    if df.empty:
        return None

    fw_metrics = residual_model_metrics(df[activity_col].values, df[fw_pred_col].values)
    ml_metrics = residual_model_metrics(df[activity_col].values, df[combined_pred_col].values)

    min_val = min(
        df[activity_col].min(),
        df[fw_pred_col].min(),
        df[combined_pred_col].min(),
    )
    max_val = max(
        df[activity_col].max(),
        df[fw_pred_col].max(),
        df[combined_pred_col].max(),
    )

    pad = (max_val - min_val) * 0.05 if max_val > min_val else 0.5
    min_val -= pad
    max_val += pad

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    sns.scatterplot(
        data=df,
        x=activity_col,
        y=fw_pred_col,
        ax=axes[0],
        s=45,
        alpha=0.75,
    )

    axes[0].plot([min_val, max_val], [min_val, max_val], "--", color="black", linewidth=1)
    axes[0].set_xlim(min_val, max_val)
    axes[0].set_ylim(min_val, max_val)
    axes[0].set_xlabel(f"Actual {activity_col}")
    axes[0].set_ylabel("FW predicted")
    axes[0].set_title(
        f"FW only\nR²={fw_metrics['r2']:.2f}, RMSE={fw_metrics['rmse']:.2f}"
    )

    sns.scatterplot(
        data=df,
        x=activity_col,
        y=combined_pred_col,
        ax=axes[1],
        s=45,
        alpha=0.75,
    )

    axes[1].plot([min_val, max_val], [min_val, max_val], "--", color="black", linewidth=1)
    axes[1].set_xlim(min_val, max_val)
    axes[1].set_ylim(min_val, max_val)
    axes[1].set_xlabel(f"Actual {activity_col}")
    axes[1].set_ylabel("FW + ML predicted")
    axes[1].set_title(
        f"FW + ML residual\nR²={ml_metrics['r2']:.2f}, RMSE={ml_metrics['rmse']:.2f}"
    )

    fig.suptitle("Observed prediction validation")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path

def mol_image_from_smiles(smiles, size=(220, 180)):
    mol = Chem.MolFromSmiles(str(smiles)) if pd.notna(smiles) else None

    if mol is None:
        return None

    try:
        return Draw.MolToImage(mol, size=size)
    except Exception:
        return None


def plot_top_candidate_cards(
    virtual_df,
    output_path,
    smiles_col="final_smiles",
    top_n=50,
    sort_col="rank_predicted_potency",
    dpi=300,
    n_cols=5,
):
    """
    Plot top virtual candidates as molecule cards with prediction metrics.
    """

    if virtual_df is None or virtual_df.empty:
        return None

    required = {
        smiles_col,
        "fw_predicted_activity",
        "ml_residual_pred",
        "combined_predicted_activity",
        "combined_pred_lower_95",
        "ml_residual_bootstrap_sd",
        "nearest_analogue_similarity",
        "applicability_domain",
        "novelty",
        "expected_improvement",
        "probability_improves_best",
        "clogP",
        "predicted_lle",
        "predicted_lle_lower_bound",
        "rank_predicted_potency",
        "rank_expected_improvement",
        "rank_confidence",
        "rank_novelty",
        "rank_lle",
        "rank_lle_confidence",
    }

    missing = required - set(virtual_df.columns)
    if missing:
        logger.warning(
            f"[FW residual ML] Candidate card plot skipped, missing columns: {missing}"
        )
        return None

    df = virtual_df.copy()

    # Keep only valid final molecules if diagnostics exist.
    if "reconstruction_success" in df.columns:
        df = df[df["reconstruction_success"] == True].copy()

    if "has_unresolved_dummy_atoms" in df.columns:
        df = df[df["has_unresolved_dummy_atoms"] == False].copy()

    numeric_card_cols = [
        "fw_predicted_activity",
        "ml_residual_pred",
        "combined_predicted_activity",
        "combined_pred_lower_95",
        "ml_residual_bootstrap_sd",
        "nearest_analogue_similarity",
        "expected_improvement",
        "probability_improves_best",
        "clogP",
        "predicted_lle",
        "predicted_lle_lower_bound",
        "rank_predicted_potency",
        "rank_expected_improvement",
        "rank_confidence",
        "rank_novelty",
        "rank_lle",
        "rank_lle_confidence",
    ]

    for col in numeric_card_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(
                df[col],
                errors="coerce",
            )

    df = df.dropna(subset=[smiles_col, sort_col])

    sort_ascending = str(sort_col).startswith("rank_")

    df = (
        df.sort_values(
            sort_col,
            ascending=sort_ascending,
        )
        .head(top_n)
        .reset_index(drop=True)
    )

    if df.empty:
        return None

    n = len(df)
    n_cols = min(n_cols, n)
    n_rows = int(np.ceil(n / n_cols))

    fig_w = max(10, n_cols * 3.0)
    fig_h = max(5, n_rows * 3.6)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_w, fig_h))
    axes = np.array(axes).reshape(-1)

    for ax in axes:
        ax.axis("off")

    for i, row in df.iterrows():
        ax = axes[i]
        ax.axis("off")

        img = mol_image_from_smiles(row[smiles_col], size=(240, 180))

        if img is not None:
            ax.imshow(img)

        title = f"Rank {i + 1}"

        if "candidate_id" in df.columns:
            title += f" | {str(row['candidate_id'])[:14]}"

        ax.set_title(title, fontsize=8)

        prob_improves = row.get(
            "probability_improves_best",
            np.nan,
        )

        prob_improves_text = (
            f"{prob_improves:.1%}"
            if pd.notna(prob_improves)
            else "NA"
        )

        txt = (
            f"FW={row['fw_predicted_activity']:.2f}\n"
            f"ML resid={row['ml_residual_pred']:+.2f}\n"
            f"Combined={row['combined_predicted_activity']:.2f}\n"
            f"95% low={row['combined_pred_lower_95']:.2f}\n"
            f"SD={row['ml_residual_bootstrap_sd']:.2f}\n"
            f"EI={row['expected_improvement']:.3f}\n"
            f"P(>best)={prob_improves_text}\n"
            f"cLogP={row['clogP']:.2f}\n"
            f"LLE={row['predicted_lle']:.2f}\n"
            f"LLE LB={row['predicted_lle_lower_bound']:.2f}\n"
            f"NN sim={row['nearest_analogue_similarity']:.2f}\n"
            f"AD={row['applicability_domain']} | "
            f"Novelty={row['novelty']}\n"
            f"Ranks: "
            f"P={int(row['rank_predicted_potency'])}, "
            f"EI={int(row['rank_expected_improvement'])}, "
            f"C={int(row['rank_confidence'])}, "
            f"N={int(row['rank_novelty'])}, "
            f"LLE={int(row['rank_lle'])}, "
            f"LLE-C={int(row['rank_lle_confidence'])}"
        )

        ax.text(
            0.02,
            -0.02,
            txt,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
        )

    fig.suptitle(f"Top {len(df)} FW + ML residual candidates ranked by {sort_col}", fontsize=14)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path

def fp_to_array(smiles, radius=2, nbits=1024):
    mol = Chem.MolFromSmiles(str(smiles)) if pd.notna(smiles) else None
    arr = np.zeros((nbits,), dtype=np.float32)

    if mol is None:
        return arr

    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def calc_physchem_descriptors(smiles):
    mol = Chem.MolFromSmiles(str(smiles)) if pd.notna(smiles) else None

    if mol is None:
        return {
            "MW": np.nan,
            "clogP": np.nan,
            "TPSA": np.nan,
            "HBD": np.nan,
            "HBA": np.nan,
            "RotBonds": np.nan,
            "Fsp3": np.nan,
            "RingCount": np.nan,
        }

    return {
        "MW": Descriptors.MolWt(mol),
        "clogP": Crippen.MolLogP(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "HBD": rdMolDescriptors.CalcNumHBD(mol),
        "HBA": rdMolDescriptors.CalcNumHBA(mol),
        "RotBonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "Fsp3": rdMolDescriptors.CalcFractionCSP3(mol),
        "RingCount": rdMolDescriptors.CalcNumRings(mol),
    }


def build_residual_ml_features(
    train_df,
    pred_df,
    r_cols,
    smiles_col_train,
    smiles_col_pred,
    fw_pred_col,
    cfg,
):
    feature_cfg = cfg.get("features", {})

    include_rgroup_onehot = feature_cfg.get("include_rgroup_onehot", True)
    include_morgan_fp = feature_cfg.get("include_morgan_fp", True)
    include_physchem = feature_cfg.get("include_physchem", True)
    include_fw_prediction = feature_cfg.get("include_fw_prediction", True)

    fp_radius = feature_cfg.get("fp_radius", 2)
    fp_nbits = feature_cfg.get("fp_nbits", 1024)

    train_parts = []
    pred_parts = []
    feature_names = []

    if include_rgroup_onehot and r_cols:
        train_r = train_df[r_cols].copy()
        pred_r = pred_df[r_cols].copy()

        for r in r_cols:
            train_r[r] = train_r[r].astype("object").where(train_r[r].notna(), "__MISSING__").astype(str)
            pred_r[r] = pred_r[r].astype("object").where(pred_r[r].notna(), "__MISSING__").astype(str)

        combined_r = pd.concat([train_r, pred_r], axis=0, ignore_index=True)
        combined_oh = pd.get_dummies(combined_r, columns=r_cols, prefix=r_cols)

        train_oh = combined_oh.iloc[:len(train_df)].reset_index(drop=True)
        pred_oh = combined_oh.iloc[len(train_df):].reset_index(drop=True)

        train_parts.append(train_oh.values.astype(np.float32))
        pred_parts.append(pred_oh.values.astype(np.float32))
        feature_names.extend(train_oh.columns.tolist())

    if include_fw_prediction:
        train_fw = pd.to_numeric(train_df[fw_pred_col], errors="coerce").fillna(train_df[fw_pred_col].median())
        pred_fw = pd.to_numeric(pred_df[fw_pred_col], errors="coerce").fillna(train_df[fw_pred_col].median())

        train_parts.append(train_fw.values.reshape(-1, 1).astype(np.float32))
        pred_parts.append(pred_fw.values.reshape(-1, 1).astype(np.float32))
        feature_names.append(fw_pred_col)

    if include_physchem:
        train_desc = pd.DataFrame([calc_physchem_descriptors(s) for s in train_df[smiles_col_train]])
        pred_desc = pd.DataFrame([calc_physchem_descriptors(s) for s in pred_df[smiles_col_pred]])

        medians = train_desc.median(numeric_only=True)
        train_desc = train_desc.fillna(medians)
        pred_desc = pred_desc.fillna(medians)

        train_parts.append(train_desc.values.astype(np.float32))
        pred_parts.append(pred_desc.values.astype(np.float32))
        feature_names.extend(train_desc.columns.tolist())

    if include_morgan_fp:
        train_fp = np.vstack([
            fp_to_array(s, radius=fp_radius, nbits=fp_nbits)
            for s in train_df[smiles_col_train]
        ])

        pred_fp = np.vstack([
            fp_to_array(s, radius=fp_radius, nbits=fp_nbits)
            for s in pred_df[smiles_col_pred]
        ])

        train_parts.append(train_fp)
        pred_parts.append(pred_fp)
        feature_names.extend([f"morgan_{i}" for i in range(fp_nbits)])

    X_train = np.hstack(train_parts) if train_parts else np.empty((len(train_df), 0))
    X_pred = np.hstack(pred_parts) if pred_parts else np.empty((len(pred_df), 0))

    return X_train, X_pred, feature_names

def make_base_residual_models(cfg):
    ensemble_cfg = cfg.get("ensemble", {})
    model_types = ensemble_cfg.get("model_types", ["random_forest", "extra_trees"])
    random_state = ensemble_cfg.get("random_state", 42)

    models = []

    if "random_forest" in model_types:
        models.append(
            RandomForestRegressor(
                n_estimators=300,
                min_samples_leaf=2,
                random_state=random_state,
                n_jobs=1,
            )
        )

    if "extra_trees" in model_types:
        models.append(
            ExtraTreesRegressor(
                n_estimators=300,
                min_samples_leaf=2,
                random_state=random_state,
                n_jobs=1,
            )
        )

    if "hist_gradient_boosting" in model_types:
        models.append(
            HistGradientBoostingRegressor(
                max_iter=300,
                learning_rate=0.05,
                random_state=random_state,
            )
        )

    return models

def _fit_single_bootstrap_model(
    bootstrap_idx,
    X,
    y,
    base_models,
    random_state,
):
    """
    Fit one independently resampled residual model.

    This function must remain at module level so that the loky
    process backend can serialize and execute it.
    """

    rng = np.random.default_rng(random_state + bootstrap_idx)

    sample_idx = rng.choice(
        len(X),
        size=len(X),
        replace=True,
    )

    base_template = base_models[bootstrap_idx % len(base_models)]
    model = clone(base_template)

    # Ensure each model also has an independent estimator seed.
    model_params = model.get_params(deep=False)

    if "random_state" in model_params:
        model.set_params(
            random_state=random_state + bootstrap_idx
        )

    model.fit(
        X[sample_idx],
        y[sample_idx],
    )

    return bootstrap_idx, model

def fit_bootstrap_residual_ensemble(X, y, cfg):
    """
    Fit bootstrap residual models in parallel.

    Parallelism is applied across bootstrap models. Individual
    RandomForest and ExtraTrees estimators use n_jobs=1 to avoid
    nested parallelism.
    """

    ensemble_cfg = cfg.get("ensemble", {})

    n_bootstrap = int(
        ensemble_cfg.get("n_bootstrap_models", 50)
    )
    random_state = int(
        ensemble_cfg.get("random_state", 42)
    )
    n_jobs = int(
        ensemble_cfg.get("n_jobs", -1)
    )
    backend = ensemble_cfg.get(
        "backend",
        "loky"
    )
    inner_max_num_threads = int(
        ensemble_cfg.get("inner_max_num_threads", 1)
    )

    if n_bootstrap < 1:
        raise ValueError(
            "[FW residual ML] n_bootstrap_models must be at least 1."
        )

    if len(X) != len(y):
        raise ValueError(
            f"[FW residual ML] X/y length mismatch: "
            f"len(X)={len(X)}, len(y)={len(y)}"
        )

    if len(X) == 0:
        raise ValueError(
            "[FW residual ML] Cannot fit bootstrap ensemble on empty data."
        )

    base_models = make_base_residual_models(cfg)

    effective_n_jobs = (
        os.cpu_count()
        if n_jobs == -1
        else max(1, min(n_jobs, n_bootstrap))
    )

    logger.debug(
        f"[FW residual ML] Starting parallel bootstrap ensemble | "
        f"models={n_bootstrap} | "
        f"workers={effective_n_jobs} | "
        f"backend={backend} | "
        f"inner_threads={inner_max_num_threads}"
    )

    with parallel_config(
        backend=backend,
        n_jobs=n_jobs,
        inner_max_num_threads=inner_max_num_threads,
    ):
        fitted = Parallel(
            verbose=0,
            batch_size=1,
        )(
            delayed(_fit_single_bootstrap_model)(
                bootstrap_idx=i,
                X=X,
                y=y,
                base_models=base_models,
                random_state=random_state,
            )
            for i in range(n_bootstrap)
        )

    # Parallel jobs may finish out of order. Restore deterministic order.
    fitted = sorted(fitted, key=lambda item: item[0])
    models = [model for _, model in fitted]

    logger.debug(
        f"[FW residual ML] Parallel bootstrap ensemble complete | "
        f"models={len(models)}"
    )

    return models

def _predict_single_residual_model(model, X):
    return model.predict(X)


def predict_residual_ensemble(models, X, cfg=None):
    """
    Predict with ensemble models in parallel across models.
    """

    if not models:
        raise ValueError(
            "[FW residual ML] Cannot predict with an empty ensemble."
        )

    if cfg is None:
        cfg = {}

    ensemble_cfg = cfg.get("ensemble", {})

    n_jobs = int(
        ensemble_cfg.get("prediction_n_jobs", ensemble_cfg.get("n_jobs", -1))
    )
    backend = ensemble_cfg.get(
        "backend",
        "loky"
    )
    inner_max_num_threads = int(
        ensemble_cfg.get("inner_max_num_threads", 1)
    )

    with parallel_config(
        backend=backend,
        n_jobs=n_jobs,
        inner_max_num_threads=inner_max_num_threads,
    ):
        pred_list = Parallel(
            verbose=0,
            batch_size=1,
        )(
            delayed(_predict_single_residual_model)(model, X)
            for model in models
        )

    preds = np.vstack(pred_list)

    return {
        "mean": preds.mean(axis=0),
        "sd": preds.std(axis=0),
        "p05": np.quantile(preds, 0.05, axis=0),
        "p95": np.quantile(preds, 0.95, axis=0),
    }

def cross_val_residual_predictions(X, y, cfg):
    validation_cfg = cfg.get("validation", {})
    cv_folds = validation_cfg.get("cv_folds", 5)
    random_state = cfg.get("ensemble", {}).get("random_state", 42)

    n = len(y)

    if n < cv_folds:
        cv_folds = max(2, n)

    kf = KFold(n_splits=cv_folds, shuffle=True, random_state=random_state)

    y_pred = np.full(n, np.nan)

    for fold_idx, (train_idx, test_idx) in enumerate(
        kf.split(X),
        start=1,
    ):
        logger.debug(
            f"[FW residual ML] "
            f"CV fold {fold_idx}/{cv_folds} | "
            f"train={len(train_idx)} | "
            f"test={len(test_idx)}"
        )
        models = fit_bootstrap_residual_ensemble(X[train_idx], y[train_idx], cfg)
        pred = predict_residual_ensemble(
            models,
            X[test_idx],
            cfg=cfg,
        )
        y_pred[test_idx] = pred["mean"]

    logger.debug(
        "[FW residual ML] Cross-validation complete"
    )

    return y_pred


def residual_model_metrics(y_true, y_pred):
    valid = ~pd.isna(y_true) & ~pd.isna(y_pred)
    y_true = np.asarray(y_true)[valid]
    y_pred = np.asarray(y_pred)[valid]

    if len(y_true) < 2:
        return {
            "n": len(y_true),
            "r2": np.nan,
            "rmse": np.nan,
            "mae": np.nan,
        }

    mse = mean_squared_error(y_true, y_pred)

    return {
        "n": len(y_true),
        "r2": r2_score(y_true, y_pred),
        "rmse": float(np.sqrt(mse)),
        "mae": mean_absolute_error(y_true, y_pred),
    }


def conformal_abs_error_quantile(y_true, y_pred, alpha=0.05):
    valid = ~pd.isna(y_true) & ~pd.isna(y_pred)
    errors = np.abs(np.asarray(y_true)[valid] - np.asarray(y_pred)[valid])

    if len(errors) == 0:
        return np.nan

    return float(np.quantile(errors, 1.0 - alpha))

def smiles_to_fp(smiles, radius=2, nbits=1024):
    mol = Chem.MolFromSmiles(str(smiles)) if pd.notna(smiles) else None

    if mol is None:
        return None

    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)


def nearest_tanimoto_similarity(query_smiles, reference_fps, radius=2, nbits=1024):
    qfp = smiles_to_fp(query_smiles, radius=radius, nbits=nbits)

    if qfp is None or not reference_fps:
        return np.nan

    sims = DataStructs.BulkTanimotoSimilarity(qfp, reference_fps)

    return float(max(sims)) if sims else np.nan


def classify_applicability_domain(similarity, bootstrap_sd, cfg):
    ad_cfg = cfg.get("applicability_domain", {})

    sim_inside = ad_cfg.get("similarity_cutoff_inside", 0.70)
    sd_cutoff = ad_cfg.get("bootstrap_sd_cutoff", 0.35)

    if pd.isna(similarity):
        return "Unknown"

    if similarity >= sim_inside and bootstrap_sd <= sd_cutoff:
        return "Inside"

    if similarity >= sim_inside:
        return "Borderline"

    return "Outside"


def classify_novelty(similarity, cfg):
    ad_cfg = cfg.get("applicability_domain", {})
    near_cutoff = ad_cfg.get("similarity_cutoff_near", 0.85)
    inside_cutoff = ad_cfg.get("similarity_cutoff_inside", 0.70)

    if pd.isna(similarity):
        return "Unknown"

    if similarity >= near_cutoff:
        return "Low"

    if similarity >= inside_cutoff:
        return "Moderate"

    return "High"

def compute_expected_improvement(
    predicted_mean,
    predicted_sd,
    best_observed,
    exploration_margin=0.0,
    min_sd=1e-8,
):
    """
    Calculate Gaussian expected improvement over the best observed activity.

    For activity values where higher is better:

        improvement = predicted_mean - best_observed - exploration_margin

        EI = improvement * Phi(z) + predicted_sd * phi(z)

        z = improvement / predicted_sd

    Parameters
    ----------
    predicted_mean : array-like
        Combined FW + ML activity prediction.

    predicted_sd : array-like
        Predictive uncertainty estimate. Here this is the bootstrap SD
        of the ML residual prediction.

    best_observed : float
        Highest observed activity within the current series/core group.

    exploration_margin : float
        Optional minimum improvement required over the current best.
        A value of 0.0 uses standard expected improvement.

    min_sd : float
        Numerical floor for SD values.

    Returns
    -------
    numpy.ndarray
        Expected improvement values.
    """

    mean = np.asarray(predicted_mean, dtype=float)
    sd = np.asarray(predicted_sd, dtype=float)

    sd_safe = np.maximum(sd, min_sd)

    improvement = (
        mean
        - float(best_observed)
        - float(exploration_margin)
    )

    z = improvement / sd_safe

    expected_improvement = (
        improvement * norm.cdf(z)
        + sd_safe * norm.pdf(z)
    )

    # Deterministic limit where SD is effectively zero.
    near_zero_sd = sd <= min_sd

    expected_improvement[near_zero_sd] = np.maximum(
        improvement[near_zero_sd],
        0.0,
    )

    expected_improvement[~np.isfinite(expected_improvement)] = np.nan

    return expected_improvement

def add_candidate_rankings(
    virtual_df,
    observed_df,
    activity_col,
    combined_pred_col="combined_predicted_activity",
    uncertainty_col="ml_residual_bootstrap_sd",
    lower_bound_col="combined_pred_lower_95",
    similarity_col="nearest_analogue_similarity",
    exploration_margin=0.0,
):
    """
    Add four independent, data-derived candidate rankings.

    Rankings
    --------
    rank_predicted_potency
        Highest combined FW + ML prediction first.

    rank_expected_improvement
        Highest expected improvement over the best observed activity first.

    rank_confidence
        Highest conservative lower prediction bound first.

    rank_novelty
        Lowest nearest-analogue similarity first. Combined predicted
        activity and lower uncertainty are used only as tie-breakers.

    No weighted priority score is calculated.
    """

    if virtual_df is None or virtual_df.empty:
        return virtual_df

    if observed_df is None or observed_df.empty:
        raise ValueError(
            "[FW residual ML] Cannot rank candidates without observed data."
        )

    candidate_df = virtual_df.copy()
    observed_work = observed_df.copy()

    required_virtual = {
        combined_pred_col,
        uncertainty_col,
        lower_bound_col,
        similarity_col,
        "predicted_lle",
        "predicted_lle_lower_bound",
    }

    missing_virtual = required_virtual - set(candidate_df.columns)

    if missing_virtual:
        raise KeyError(
            f"[FW residual ML] Cannot calculate candidate rankings. "
            f"Missing virtual columns: {sorted(missing_virtual)}"
        )

    if activity_col not in observed_work.columns:
        raise KeyError(
            f"[FW residual ML] Observed activity column "
            f"'{activity_col}' is missing."
        )

    numeric_virtual_cols = [
        combined_pred_col,
        uncertainty_col,
        lower_bound_col,
        similarity_col,
        "predicted_lle",
        "predicted_lle_lower_bound",
    ]

    for col in numeric_virtual_cols:
        candidate_df[col] = pd.to_numeric(
            candidate_df[col],
            errors="coerce",
        )

    observed_activity = pd.to_numeric(
        observed_work[activity_col],
        errors="coerce",
    ).dropna()

    if observed_activity.empty:
        raise ValueError(
            "[FW residual ML] No valid observed activity values "
            "available for candidate ranking."
        )

    best_observed = float(observed_activity.max())

    candidate_df["best_observed_activity"] = best_observed

    candidate_df["predicted_improvement_over_best"] = (
        candidate_df[combined_pred_col]
        - best_observed
    )

    candidate_df["expected_improvement"] = compute_expected_improvement(
        predicted_mean=candidate_df[combined_pred_col].values,
        predicted_sd=candidate_df[uncertainty_col].values,
        best_observed=best_observed,
        exploration_margin=exploration_margin,
    )

    candidate_df["probability_improves_best"] = np.nan

    valid_uncertainty = (
        candidate_df[uncertainty_col].notna()
        & (candidate_df[uncertainty_col] > 0)
        & candidate_df[combined_pred_col].notna()
    )

    candidate_df.loc[
        valid_uncertainty,
        "probability_improves_best",
    ] = norm.cdf(
        (
            candidate_df.loc[
                valid_uncertainty,
                combined_pred_col,
            ]
            - best_observed
            - exploration_margin
        )
        / candidate_df.loc[
            valid_uncertainty,
            uncertainty_col,
        ]
    )

    deterministic_rows = (
        candidate_df[uncertainty_col].notna()
        & (candidate_df[uncertainty_col] <= 0)
        & candidate_df[combined_pred_col].notna()
    )

    candidate_df.loc[
        deterministic_rows,
        "probability_improves_best",
    ] = (
        candidate_df.loc[
            deterministic_rows,
            combined_pred_col,
        ]
        > best_observed + exploration_margin
    ).astype(float)

    # Rank 1: highest combined predicted potency.
    candidate_df["rank_predicted_potency"] = (
        candidate_df[combined_pred_col]
        .rank(
            method="min",
            ascending=False,
            na_option="bottom",
        )
        .astype("Int64")
    )

    # Rank 2: highest expected improvement.
    candidate_df["rank_expected_improvement"] = (
        candidate_df["expected_improvement"]
        .rank(
            method="min",
            ascending=False,
            na_option="bottom",
        )
        .astype("Int64")
    )

    # Rank 3: highest conservative lower bound.
    candidate_df["rank_confidence"] = (
        candidate_df[lower_bound_col]
        .rank(
            method="min",
            ascending=False,
            na_option="bottom",
        )
        .astype("Int64")
    )

    # Rank 4: highest conservative lower bound.

    candidate_df["rank_lle"] = (
        candidate_df["predicted_lle"]
        .rank(
            method="min",
            ascending=False,
            na_option="bottom",
        )
        .astype("Int64")
    )

    candidate_df["rank_lle_confidence"] = (
        candidate_df["predicted_lle_lower_bound"]
        .rank(
            method="min",
            ascending=False,
            na_option="bottom",
        )
        .astype("Int64")
    )

    # Rank 5: most novel first.
    #
    # This uses lexicographic ordering rather than a weighted score:
    #   1. lower nearest-neighbour similarity,
    #   2. higher combined predicted activity,
    #   3. lower bootstrap uncertainty.
    #
    # Therefore novelty is never converted into an arbitrary activity bonus.
    novelty_order = candidate_df.sort_values(
        by=[
            similarity_col,
            combined_pred_col,
            uncertainty_col,
        ],
        ascending=[
            True,
            False,
            True,
        ],
        na_position="last",
        kind="mergesort",
    ).index

    novelty_rank = pd.Series(
        np.arange(1, len(candidate_df) + 1),
        index=novelty_order,
        dtype="int64",
    )

    candidate_df["rank_novelty"] = (
        novelty_rank
        .reindex(candidate_df.index)
        .astype("Int64")
    )

    n_candidates = len(candidate_df)

    for rank_col in [
        "rank_predicted_potency",
        "rank_expected_improvement",
        "rank_confidence",
        "rank_novelty",
        "rank_lle",
        "rank_lle_confidence",
    ]:
        if rank_col not in candidate_df.columns:
            continue

        pct_col = rank_col.replace("rank_", "pct_")

        candidate_df[pct_col] = (
            100.0
            * (
                1.0
                - (
                    candidate_df[rank_col] - 1
                )
                / max(1, n_candidates - 1)
            )
        )


    return candidate_df

def run_fw_residual_ml_for_group(
    observed_df,
    virtual_df,
    cfg,
    output_dir,
    series,
    core_hash,
):
    output_dir.mkdir(parents=True, exist_ok=True)

    activity_col = cfg.get("activity_col", "pIC50")
    fw_pred_col = cfg.get("fw_pred_col", "fw_predicted_activity")
    residual_col = cfg.get("residual_col", "fw_residual")
    smiles_col_observed = cfg.get("smiles_col_observed", "canonical_input_smiles")
    smiles_col_virtual = cfg.get("smiles_col_virtual", "final_smiles")

    min_train = cfg.get("min_train_compounds_per_group", 25)
    min_virtual = cfg.get("min_virtual_candidates_per_group", 1)

    observed_df = observed_df.copy()
    virtual_df = virtual_df.copy()

    if len(observed_df) < min_train:
        logger.info(
            f"[FW residual ML] Skipping series={series}, core={core_hash}: "
            f"n_train={len(observed_df)} < {min_train}"
        )
        return None

    if len(virtual_df) < min_virtual:
        logger.info(
            f"[FW residual ML] Skipping series={series}, core={core_hash}: "
            f"n_virtual={len(virtual_df)} < {min_virtual}"
        )
        return None

    # Keep only valid virtual final molecules.
    if "reconstruction_success" in virtual_df.columns:
        virtual_df = virtual_df[virtual_df["reconstruction_success"] == True].copy()

    if "has_unresolved_dummy_atoms" in virtual_df.columns:
        virtual_df = virtual_df[virtual_df["has_unresolved_dummy_atoms"] == False].copy()

    virtual_df = virtual_df.dropna(subset=[smiles_col_virtual, fw_pred_col])

    if virtual_df.empty:
        logger.info(
            f"[FW residual ML] Skipping series={series}, core={core_hash}: "
            f"no valid virtual molecules."
        )
        return None

    observed_df[activity_col] = pd.to_numeric(observed_df[activity_col], errors="coerce")
    observed_df[fw_pred_col] = pd.to_numeric(observed_df[fw_pred_col], errors="coerce")

    if residual_col not in observed_df.columns:
        observed_df[residual_col] = observed_df[activity_col] - observed_df[fw_pred_col]

    observed_df[residual_col] = pd.to_numeric(observed_df[residual_col], errors="coerce")

    observed_df = observed_df.dropna(
        subset=[activity_col, fw_pred_col, residual_col, smiles_col_observed]
    ).copy()

    if len(observed_df) < min_train:
        return None

    r_cols = get_rgroup_columns(observed_df)
    r_cols = [r for r in r_cols if r in virtual_df.columns]

    logger.info(
        f"[FW residual ML] Building features "
        f"for series={series}, core={core_hash} | "
        f"observed={len(observed_df)} | "
        f"virtual={len(virtual_df)}"
    )

    X_train, X_virtual, feature_names = build_residual_ml_features(
        train_df=observed_df,
        pred_df=virtual_df,
        r_cols=r_cols,
        smiles_col_train=smiles_col_observed,
        smiles_col_pred=smiles_col_virtual,
        fw_pred_col=fw_pred_col,
        cfg=cfg,
    )

    logger.info(
        f"[FW residual ML] Feature matrix built "
        f"for series={series}, core={core_hash} | "
        f"train_shape={X_train.shape} | "
        f"virtual_shape={X_virtual.shape}"
    )

    y_resid = observed_df[residual_col].values.astype(float)

    cv_pred_resid = cross_val_residual_predictions(X_train, y_resid, cfg)

    logger.info(
        f"[FW residual ML] CV residual predictions complete "
        f"for series={series}, core={core_hash}"
    )

    observed_df["ml_residual_pred_cv"] = cv_pred_resid
    observed_df["combined_predicted_activity_cv"] = (
        observed_df[fw_pred_col] + observed_df["ml_residual_pred_cv"]
    )

    fw_metrics = residual_model_metrics(
        observed_df[activity_col].values,
        observed_df[fw_pred_col].values,
    )

    combined_metrics = residual_model_metrics(
        observed_df[activity_col].values,
        observed_df["combined_predicted_activity_cv"].values,
    )

    delta_rmse = (
        fw_metrics["rmse"] - combined_metrics["rmse"]
    )

    delta_r2 = (
        combined_metrics["r2"] - fw_metrics["r2"]
    )

    logger.info(
        f"[FW residual ML] "
        f"series={series}, core={core_hash} | "
        f"FW R2={fw_metrics['r2']:.3f} | "
        f"FW+ML R2={combined_metrics['r2']:.3f} | "
        f"ΔR2={delta_r2:+.3f} | "
        f"ΔRMSE={delta_rmse:+.3f}"
    )

    residual_metrics = residual_model_metrics(
        y_resid,
        cv_pred_resid,
    )

    alpha = cfg.get("validation", {}).get("interval_alpha", 0.05)
    conformal_q = conformal_abs_error_quantile(y_resid, cv_pred_resid, alpha=alpha)

    logger.info(
            f"[FW residual ML] Training final ensemble "
            f"for series={series}, core={core_hash} | "
            f"n_train={len(X_train)}"
        )

    models = fit_bootstrap_residual_ensemble(X_train, y_resid, cfg)

    logger.info(
        f"[FW residual ML] Final ensemble trained "
        f"for series={series}, core={core_hash}"
    )

    obs_ens = predict_residual_ensemble(
        models,
        X_train,
        cfg=cfg,
    )

    logger.debug(
        f"[FW residual ML] Scoring "
        f"{len(virtual_df):,} virtual molecules "
        f"for series={series}, core={core_hash}"
    )

    virt_ens = predict_residual_ensemble(
        models,
        X_virtual,
        cfg=cfg,
    )

    logger.info(
        f"[FW residual ML] Virtual scoring complete "
        f"for series={series}, core={core_hash}"
    )

    observed_df["ml_residual_pred"] = obs_ens["mean"]
    observed_df["ml_residual_bootstrap_sd"] = obs_ens["sd"]
    observed_df["combined_predicted_activity"] = (
        observed_df[fw_pred_col] + observed_df["ml_residual_pred"]
    )

    virtual_df["ml_residual_pred"] = virt_ens["mean"]
    virtual_df["ml_residual_bootstrap_sd"] = virt_ens["sd"]
    virtual_df["ml_residual_p05"] = virt_ens["p05"]
    virtual_df["ml_residual_p95"] = virt_ens["p95"]

    virtual_df["combined_predicted_activity"] = (
        virtual_df[fw_pred_col]
        + virtual_df["ml_residual_pred"]
    )

    # Add conformal interval before calculating lower-bound LLE.
    virtual_df["conformal_abs_error_q95"] = conformal_q

    virtual_df["combined_pred_lower_95"] = (
        virtual_df["combined_predicted_activity"]
        - conformal_q
    )

    virtual_df["combined_pred_upper_95"] = (
        virtual_df["combined_predicted_activity"]
        + conformal_q
    )

    # Calculate and retain physicochemical descriptors.
    descriptor_df = pd.DataFrame(
        [
            calc_physchem_descriptors(smiles)
            for smiles in virtual_df[smiles_col_virtual]
        ],
        index=virtual_df.index,
    )

    for col in descriptor_df.columns:
        virtual_df[col] = descriptor_df[col]

    # Ensure the descriptor needed for LLE is numeric.
    virtual_df["clogP"] = pd.to_numeric(
        virtual_df["clogP"],
        errors="coerce",
    )

    # Predicted LLE based on the mean combined potency prediction.
    virtual_df["predicted_lle"] = (
        virtual_df["combined_predicted_activity"]
        - virtual_df["clogP"]
    )

    # Conservative LLE based on the lower potency bound.
    virtual_df["predicted_lle_lower_bound"] = (
        virtual_df["combined_pred_lower_95"]
        - virtual_df["clogP"]
    )

    descriptor_df = pd.DataFrame([
        calc_physchem_descriptors(s)
        for s in virtual_df[smiles_col_virtual]
    ])

    for col in descriptor_df.columns:
        virtual_df[col] = descriptor_df[col].values

    virtual_df["conformal_abs_error_q95"] = conformal_q
    virtual_df["combined_pred_lower_95"] = (
        virtual_df["combined_predicted_activity"] - conformal_q
    )
    virtual_df["combined_pred_upper_95"] = (
        virtual_df["combined_predicted_activity"] + conformal_q
    )

    fp_radius = cfg.get("features", {}).get("fp_radius", 2)
    fp_nbits = cfg.get("features", {}).get("fp_nbits", 1024)

    reference_fps = [
        smiles_to_fp(s, radius=fp_radius, nbits=fp_nbits)
        for s in observed_df[smiles_col_observed]
    ]
    reference_fps = [fp for fp in reference_fps if fp is not None]

    virtual_df["nearest_analogue_similarity"] = virtual_df[smiles_col_virtual].apply(
        lambda s: nearest_tanimoto_similarity(
            s,
            reference_fps,
            radius=fp_radius,
            nbits=fp_nbits,
        )
    )

    virtual_df["applicability_domain"] = virtual_df.apply(
        lambda row: classify_applicability_domain(
            row["nearest_analogue_similarity"],
            row["ml_residual_bootstrap_sd"],
            cfg,
        ),
        axis=1,
    )

    virtual_df["novelty"] = virtual_df["nearest_analogue_similarity"].apply(
        lambda sim: classify_novelty(sim, cfg)
    )

    ranking_cfg = cfg.get("ranking", {})

    exploration_margin = ranking_cfg.get(
        "expected_improvement_margin",
        0.0,
    )

    virtual_df = add_candidate_rankings(
        virtual_df=virtual_df,
        observed_df=observed_df,
        activity_col=activity_col,
        combined_pred_col="combined_predicted_activity",
        uncertainty_col="ml_residual_bootstrap_sd",
        lower_bound_col="combined_pred_lower_95",
        similarity_col="nearest_analogue_similarity",
        exploration_margin=exploration_margin,
    )

    # Remove the legacy score if this task is rerun using an input dataframe
    # that already contains it.
    virtual_df = virtual_df.drop(
        columns=["priority_score"],
        errors="ignore",
    )

    top_n = ranking_cfg.get("top_n_per_group", 100)

    # Keep a union of the top candidates from every ranking.
    # A candidate appearing in several rankings is retained once.
    top_candidate_indices = set()

    for rank_col in [
        "rank_predicted_potency",
        "rank_expected_improvement",
        "rank_confidence",
        "rank_novelty",
    ]:
        selected_indices = (
            virtual_df
            .nsmallest(top_n, rank_col)
            .index
            .tolist()
        )

        top_candidate_indices.update(selected_indices)

    top_df = virtual_df.loc[
        sorted(top_candidate_indices)
    ].copy()

    # Sort the combined review table by potency rank only for stable display.
    # The other independent ranks remain available as columns.
    top_df = top_df.sort_values(
        [
            "rank_predicted_potency",
            "rank_expected_improvement",
            "rank_confidence",
            "rank_lle",
            "rank_lle_confidence",
            "rank_novelty",
        ],
        ascending=True,
    ).reset_index(drop=True)

    observed_csv = output_dir / "observed_residual_predictions.csv"
    virtual_csv = output_dir / "virtual_residual_predictions.csv"
    top_csv = output_dir / "top_candidates.csv"

    observed_df.to_csv(observed_csv, index=False)
    virtual_df.to_csv(virtual_csv, index=False)
    top_df.to_csv(top_csv, index=False)
    
    diagnostics_cfg = cfg.get("diagnostics", {})
    make_plots = diagnostics_cfg.get("make_plots", True)

    observed_validation_plot = None
    virtual_components_plot = None
    candidate_card_plots = {}

    if make_plots:
        observed_validation_plot = output_dir / "observed_fw_vs_ml_validation.png"
        plot_observed_fw_vs_ml_validation(
            observed_df=observed_df,
            output_path=observed_validation_plot,
            activity_col=activity_col,
            fw_pred_col=fw_pred_col,
            combined_pred_col="combined_predicted_activity_cv",
        )

        virtual_components_plot = output_dir / "virtual_prediction_components.png"
        plot_virtual_prediction_components(
            virtual_df=virtual_df,
            output_path=virtual_components_plot,
            top_n=diagnostics_cfg.get("component_plot_top_n", 100),
            sort_col=diagnostics_cfg.get(
                "component_plot_sort_col",
                "rank_predicted_potency",
            ),
        )

        ranking_plot_specs = {
            "potency": "rank_predicted_potency",
            "expected_improvement": "rank_expected_improvement",
            "confidence": "rank_confidence",
            "novelty": "rank_novelty",
            "lle": "rank_lle",
            "lle_confidence": "rank_lle_confidence",
        }

        for ranking_name, rank_col in ranking_plot_specs.items():
            card_path = (
                output_dir
                / f"top_candidate_cards_{ranking_name}.png"
            )

            made_plot = plot_top_candidate_cards(
                virtual_df=virtual_df,
                output_path=card_path,
                smiles_col=smiles_col_virtual,
                top_n=diagnostics_cfg.get(
                    "candidate_card_top_n",
                    50,
                ),
                sort_col=rank_col,
                n_cols=diagnostics_cfg.get(
                    "candidate_card_n_cols",
                    5,
                ),
            )

            candidate_card_plots[ranking_name] = (
                str(made_plot)
                if made_plot is not None
                else None
            )

    summary = {
        "series": series,
        "core_hash": core_hash,
        "n_train": len(observed_df),
        "n_virtual_valid": len(virtual_df),
        "n_features": X_train.shape[1],
        "conformal_abs_error_q95": conformal_q,
        "fw_r2": fw_metrics["r2"],
        "fw_rmse": fw_metrics["rmse"],
        "fw_mae": fw_metrics["mae"],
        "combined_cv_r2": combined_metrics["r2"],
        "combined_cv_rmse": combined_metrics["rmse"],
        "combined_cv_mae": combined_metrics["mae"],
        "residual_cv_r2": residual_metrics["r2"],
        "residual_cv_rmse": residual_metrics["rmse"],
        "residual_cv_mae": residual_metrics["mae"],
        "observed_csv": str(observed_csv),
        "virtual_csv": str(virtual_csv),
        "top_candidates_csv": str(top_csv),
        "observed_validation_plot": str(observed_validation_plot) if observed_validation_plot else None,
        "virtual_components_plot": str(virtual_components_plot) if virtual_components_plot else None,
        "top_candidate_cards_potency": candidate_card_plots.get("potency"),
        "top_candidate_cards_expected_improvement": candidate_card_plots.get(
            "expected_improvement"
        ),
        "top_candidate_cards_confidence": candidate_card_plots.get("confidence"),
        "top_candidate_cards_novelty": candidate_card_plots.get("novelty"),
        "top_candidate_cards_lle": candidate_card_plots.get("lle"),
        "top_candidate_cards_lle_confidence": candidate_card_plots.get(
            "lle_confidence"
        ),
        "best_observed_activity": (
            float(virtual_df["best_observed_activity"].iloc[0])
            if (
                not virtual_df.empty
                and "best_observed_activity" in virtual_df.columns
            )
            else np.nan
        ),
        "max_expected_improvement": (
            float(virtual_df["expected_improvement"].max())
            if (
                not virtual_df.empty
                and "expected_improvement" in virtual_df.columns
            )
            else np.nan
        ),
    }

    logger.info(
        f"[FW residual ML] COMPLETE | "
        f"series={series} | "
        f"core={core_hash} | "
        f"n_train={len(observed_df)} | "
        f"n_virtual={len(virtual_df)} | "
        f"FW_R2={fw_metrics['r2']:.2f} | "
        f"FWML_R2={combined_metrics['r2']:.2f}"
    )

    return {
        "observed": observed_df,
        "virtual": virtual_df,
        "top": top_df,
        "summary": summary,
    }