import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw, Descriptors, Crippen, rdMolDescriptors, QED, rdRGroupDecomposition, Scaffolds

from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor, HistGradientBoostingRegressor
from sklearn.model_selection import KFold
from sklearn.base import clone
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from rdkit.Chem import AllChem, DataStructs, Descriptors, Crippen, rdMolDescriptors

from modules.utils.free_wilson_helpers import get_rgroup_columns

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

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
                n_jobs=-1,
            )
        )

    if "extra_trees" in model_types:
        models.append(
            ExtraTreesRegressor(
                n_estimators=300,
                min_samples_leaf=2,
                random_state=random_state,
                n_jobs=-1,
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

def fit_bootstrap_residual_ensemble(X, y, cfg):
    ensemble_cfg = cfg.get("ensemble", {})
    n_bootstrap = ensemble_cfg.get("n_bootstrap_models", 50)
    random_state = ensemble_cfg.get("random_state", 42)

    rng = np.random.default_rng(random_state)
    base_models = make_base_residual_models(cfg)

    models = []

    for i in range(n_bootstrap):
        if (
            i == 0
            or (i + 1) % 10 == 0
            or (i + 1) == n_bootstrap
        ):
            logger.info(
                f"[FW residual ML] Bootstrap "
                f"{i + 1}/{n_bootstrap}"
            )
        idx = rng.choice(len(X), size=len(X), replace=True)
        base = clone(base_models[i % len(base_models)])
        base.fit(X[idx], y[idx])
        models.append(base)

    return models


def predict_residual_ensemble(models, X):
    preds = np.vstack([m.predict(X) for m in models])

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
        logger.info(
            f"[FW residual ML] "
            f"CV fold {fold_idx}/{cv_folds} | "
            f"train={len(train_idx)} | "
            f"test={len(test_idx)}"
        )
        models = fit_bootstrap_residual_ensemble(X[train_idx], y[train_idx], cfg)
        pred = predict_residual_ensemble(models, X[test_idx])
        y_pred[test_idx] = pred["mean"]

    logger.info(
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

    models = fit_bootstrap_residual_ensemble(X_train, y_resid, cfg)

    logger.info(
        f"[FW residual ML] Final ensemble trained "
        f"for series={series}, core={core_hash}"
    )

    obs_ens = predict_residual_ensemble(models, X_train)

    logger.info(
        f"[FW residual ML] Scoring "
        f"{len(virtual_df):,} virtual molecules "
        f"for series={series}, core={core_hash}"
    )

    virt_ens = predict_residual_ensemble(models, X_virtual)

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
        virtual_df[fw_pred_col] + virtual_df["ml_residual_pred"]
    )

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
    uncertainty_weight = ranking_cfg.get("uncertainty_weight", 0.5)
    outside_penalty = ranking_cfg.get("outside_ad_penalty", 0.5)
    novelty_bonus_moderate = ranking_cfg.get("novelty_bonus_moderate", 0.1)
    novelty_penalty_high = ranking_cfg.get("novelty_penalty_high", 0.1)

    virtual_df["priority_score"] = virtual_df["combined_predicted_activity"]
    virtual_df["priority_score"] -= uncertainty_weight * virtual_df["ml_residual_bootstrap_sd"]

    virtual_df.loc[
        virtual_df["applicability_domain"] == "Outside",
        "priority_score"
    ] -= outside_penalty

    virtual_df.loc[
        virtual_df["novelty"] == "Moderate",
        "priority_score"
    ] += novelty_bonus_moderate

    virtual_df.loc[
        virtual_df["novelty"] == "High",
        "priority_score"
    ] -= novelty_penalty_high

    virtual_df = virtual_df.sort_values("priority_score", ascending=False).reset_index(drop=True)

    top_n = ranking_cfg.get("top_n_per_group", 100)
    top_df = virtual_df.head(top_n).copy()

    observed_csv = output_dir / "observed_residual_predictions.csv"
    virtual_csv = output_dir / "virtual_residual_predictions.csv"
    top_csv = output_dir / "top_candidates.csv"

    observed_df.to_csv(observed_csv, index=False)
    virtual_df.to_csv(virtual_csv, index=False)
    top_df.to_csv(top_csv, index=False)

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