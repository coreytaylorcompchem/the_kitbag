from pipeline.task_registry import register_task

import os
import inspect
import json

from pathlib import Path
import pandas as pd
import numpy as np
from tqdm import tqdm

from scipy.stats import sem
from scipy.integrate import solve_ivp
import yaml

import matplotlib.pyplot as plt
import seaborn as sns
from cycler import cycler

import torch
from torch_geometric.loader import DataLoader

from rdkit import Chem
from rdkit.Chem import Descriptors
import joblib

from pipeline.task_registry import register_task

#Load ML models

from models.clint.model_def import ClintModel
from models.logbb.model_def import LogBBModel
from models.ppb_fu.model_def import PPBFuModel

from models.clint.featurisation import mol_to_graph as mol_to_graph_clint
from models.logbb.featurisation import mol_to_graph as mol_to_graph_logbb
from models.ppb_fu.featurisation import mol_to_graph as mol_to_graph_ppbfu

from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

def get_available_target_models():
    """
    Returns mapping of target-specific efficacy models to their
    checkpoints, classes, and featurisers.
    """
    model_dir = Path(__file__).resolve().parent.parent / "models" / "targets"

    models = {}
    for target_dir in model_dir.glob("*"):
        model_file = target_dir / "activity_gin.pt"
        if model_file.exists():
            module_name = f"models.targets.{target_dir.name}.model_def"
            featuriser_module = f"models.targets.{target_dir.name}.featurisation"

            model_class = getattr(__import__(module_name, fromlist=["ActivityModel"]), "ActivityModel")
            featuriser = getattr(__import__(featuriser_module, fromlist=["mol_to_graph"]), "mol_to_graph")

            models[target_dir.name] = {
                "path": model_file,
                "class": model_class,
                "featuriser": featuriser,
            }

    return models

@register_task("pd_parameter_prediction", category="PBPK-PD", description="Predict target-level EC50 or Emax parameters")
def pd_parameter_prediction(config, df=None):
    pd_cfg = config.get("pd_parameter_prediction", {})
    target_cfg = pd_cfg.get("targets", {})  # e.g., {"AGTR2": True}
    output_dir = Path(pd_cfg.get("output_dir", "outputs/pd/params"))
    output_dir.mkdir(parents=True, exist_ok=True)

    all_targets = get_available_target_models()
    selected = {t: all_targets[t] for t, enabled in target_cfg.items() if enabled and t in all_targets}

    if not selected:
        raise ValueError("No PD targets enabled in YAML config.")

    df = df or pd.read_csv(config["input_file"])
    df["smiles"] = df["smiles"].astype(str).str.strip()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    for target, info in selected.items():
        model = info["class"]()
        checkpoint = torch.load(info["path"], map_location=device)
        model.load_state_dict(checkpoint["model_state_dict"])
        model.to(device)
        model.eval()

        featuriser = info["featuriser"]
        preds = []
        for smi in df["smiles"]:
            graph = featuriser(smi)
            if graph is None:
                preds.append(None)
                continue
            with torch.no_grad():
                out = model(graph.to(device))
                preds.append(out.item())
        df[f"EC50_{target}"] = preds

    df.to_csv(output_dir / "pd_params.csv", index=False)
    return {"df": df}

@register_task("pbpk_pd_simulation", category="PBPK-PD",
               description="Simulate receptor inhibition using predicted IC50 and PK profiles.")
def pbpk_pd_simulation(config, df=None, **kwargs):
    """
    Computes receptor inhibition E(t) over time using the inhibitory Emax model.
    Predicts doses required to reach specified inhibition targets over given dosing intervals.
    Generates visualizations for % inhibition curves and dose vs. inhibition heatmaps.
    """
    if df is None:
        raise ValueError("Missing 'df' containing PK results + IC50 predictions.")

    pd_cfg = config.get("pbpk_pd_simulation", {})

    # Basic PD parameters
    E0 = pd_cfg.get("E0", 1.0)
    Imax = pd_cfg.get("Imax", 1.0)
    gamma = pd_cfg.get("gamma", 1.0)

    # Dosing / target info
    target_inhibitions = pd_cfg.get("target_inhibitions", [95])
    dosing_intervals = pd_cfg.get("dosing_intervals", [12])
    
    # Fields in df
    concentration_field = pd_cfg.get("concentration_field", "Cmax_plasma")
    fu_field = pd_cfg.get("fu_field", "fu_p")

    # Output options
    output_dir = Path(pd_cfg.get("output_dir", "outputs/pbpk/pd"))
    output_dir.mkdir(parents=True, exist_ok=True)
    save_heatmaps = pd_cfg.get("save_heatmaps", True)
    save_curves = pd_cfg.get("save_curves", True)

    # Dose scanning
    dose_scan_cfg = pd_cfg.get("dose_scan", {})
    dose_scan_enabled = dose_scan_cfg.get("enabled", True)
    dose_step = dose_scan_cfg.get("step_uM", 0.1)
    dose_max = dose_scan_cfg.get("max_uM", 50)

    all_results = []

    for _, row in df.iterrows():
        compound_name = row.get("compound", row.get("name", "unknown"))
        IC50 = row.get("pred_IC50_AGTR2_uM", np.nan)
        fu_p = max(row.get(fu_field, 0.1), 1e-3)
        Cmax_p = row.get(concentration_field, np.nan)

        if np.isnan(IC50) or np.isnan(Cmax_p):
            continue

        # Free concentration
        C_free = fu_p * Cmax_p

        # Compute basic inhibition
        E = E0 - Imax * (C_free ** gamma) / (IC50 ** gamma + C_free ** gamma)
        inhibition_pct = 100 * (1 - E / E0)

        # Prepare results per compound
        compound_results = {
            "compound": compound_name,
            "IC50_uM": IC50,
            "Cmax_free_uM": C_free,
            "predicted_inhibition_%": inhibition_pct,
            "dose_predictions": {},
        }

        # Dose scanning to reach target inhibitions over dosing intervals
        if dose_scan_enabled:
            for target in target_inhibitions:
                target_frac = 1 - target / 100
                predicted_doses = {}
                for interval in dosing_intervals:
                    # Using simple Emax formula solved for concentration
                    try:
                        C_required = IC50 * ((Imax / (I0 := E0*(1-target_frac)) - 1) ** (1 / gamma))
                        # Scale dose proportionally (assumes linear PK for simplicity)
                        predicted_dose = min(max(C_required / fu_p, 0), dose_max)
                    except Exception:
                        predicted_dose = np.nan
                    predicted_doses[interval] = predicted_dose
                compound_results["dose_predictions"][target] = predicted_doses

        all_results.append(compound_results)

        # Plot % inhibition curve if desired
        if save_curves:
            time_h = np.linspace(0, max(dosing_intervals), 100)
            E_t = E0 - Imax * (C_free ** gamma) / (IC50 ** gamma + C_free ** gamma)
            plt.figure()
            plt.plot(time_h, 100*(1-E_t/E0))
            plt.title(f"% Inhibition vs Time: {compound_name}")
            plt.xlabel("Time (h)")
            plt.ylabel("% Inhibition")
            plt.ylim(0, 100)
            plt.grid(True)
            plt.savefig(output_dir / f"{compound_name}_inhibition_curve.png")
            plt.close()

    # Save a summary CSV
    summary_rows = []
    for res in all_results:
        for target, doses in res["dose_predictions"].items():
            for interval, dose in doses.items():
                summary_rows.append({
                    "compound": res["compound"],
                    "IC50_uM": res["IC50_uM"],
                    "Cmax_free_uM": res["Cmax_free_uM"],
                    "target_inhibition_%": target,
                    "dosing_interval_h": interval,
                    "predicted_dose_uM": dose
                })
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(output_dir / "pd_dose_predictions.csv", index=False)
    logger.info(f"PD simulation & dose predictions saved to {output_dir / 'pd_dose_predictions.csv'}")

    # Optional heatmap
    if save_heatmaps and not summary_df.empty:
        plt.figure(figsize=(10, 6))
        pivot_table = summary_df.pivot("compound", "dosing_interval_h", "predicted_dose_uM")
        sns.heatmap(pivot_table, annot=True, fmt=".2f", cmap="YlGnBu")
        plt.title("Predicted Dose (uM) vs Compound & Dosing Interval")
        plt.ylabel("Compound")
        plt.xlabel("Dosing Interval (h)")
        plt.tight_layout()
        plt.savefig(output_dir / "dose_prediction_heatmap.png")
        plt.close()

    return {"df": summary_df, "full_results": all_results}


@register_task("pkpd_analysis", category="PBPK-PD",
               description="Analyse PK–PD simulation outputs to compute exposure–response metrics.")
def pkpd_analysis(config, df=None, **kwargs):
    """
    Post-process PD simulation outputs (from pbpk_pd_simulation).
    Computes key exposure–response metrics such as:
        - Emax (maximum inhibition)
        - AUEC (area under effect–time curve)
        - T>50% (time above 50% inhibition)
        - Mean ± SD of inhibition across virtual subjects (if population data)
    """
    if df is None:
        raise ValueError("Missing PD simulation dataframe for analysis.")

    analysis_cfg = config.get("pkpd_analysis", {})
    output_dir = Path(analysis_cfg.get("output_dir", "outputs/pbpk/pkpd_analysis"))
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Starting PK–PD analysis. Using config: {analysis_cfg}")

    # ----------------------------------------------------------------------
    # Handle both static and time-dependent results
    # ----------------------------------------------------------------------
    time_col = "time" if "time" in df.columns else None
    target_col = "target" if "target" in df.columns else None
    subject_col = "subject" if "subject" in df.columns else None

    # Convert inhibition fraction if needed
    if "inhibition_%" in df.columns:
        df["inhibition_frac"] = df["inhibition_%"] / 100.0
    elif "predicted_inhibition_%" in df.columns:
        df["inhibition_frac"] = df["predicted_inhibition_%"] / 100.0
    else:
        raise ValueError("No inhibition column found in PD dataframe.")

    results = []

    # ----------------------------------------------------------------------
    # Time-dependent aggregation (if Cp(t) → inhibition(t))
    # ----------------------------------------------------------------------
    if time_col:
        for (compound, group_df) in df.groupby("compound"):
            target = group_df[target_col].iloc[0] if target_col else None

            # Aggregate per-subject if present
            if subject_col:
                subj_metrics = []
                for subject, subj_df in group_df.groupby(subject_col):
                    times = subj_df[time_col].values
                    inhib = subj_df["inhibition_frac"].values

                    emax = inhib.max()
                    auec = np.trapz(inhib, times)
                    t50 = np.trapz(inhib >= 0.5, times)  # time above 50% inhibition

                    subj_metrics.append((emax, auec, t50))

                subj_metrics = np.array(subj_metrics)
                results.append({
                    "compound": compound,
                    "target": target,
                    "Emax_mean": subj_metrics[:,0].mean(),
                    "Emax_sd": subj_metrics[:,0].std(),
                    "AUEC_mean": subj_metrics[:,1].mean(),
                    "AUEC_sd": subj_metrics[:,1].std(),
                    "T>50%_mean": subj_metrics[:,2].mean(),
                    "T>50%_sd": subj_metrics[:,2].std(),
                })
            else:
                # Single time-course (no population)
                times = group_df[time_col].values
                inhib = group_df["inhibition_frac"].values
                emax = inhib.max()
                auec = np.trapz(inhib, times)
                t50 = np.trapz(inhib >= 0.5, times)
                results.append({
                    "compound": compound,
                    "target": target,
                    "Emax": emax,
                    "AUEC": auec,
                    "T>50%": t50,
                })
    else:
        # ----------------------------------------------------------------------
        # Static results (single inhibition value per compound)
        # ----------------------------------------------------------------------
        for _, row in df.iterrows():
            compound = row.get("compound", "unknown")
            target = row.get(target_col, None)
            results.append({
                "compound": compound,
                "target": target,
                "predicted_inhibition_%": row.get("predicted_inhibition_%", np.nan),
                "IC50_uM": row.get("IC50_uM", np.nan),
                "Cmax_free_uM": row.get("Cmax_free_uM", np.nan)
            })

    result_df = pd.DataFrame(results)

    # Save results
    summary_path = output_dir / "pkpd_summary.csv"
    result_df.to_csv(summary_path, index=False)
    logger.info(f"PK–PD analysis summary saved to {summary_path}")

    # ----------------------------------------------------------------------
    # Plotting inhibition curves
    # ----------------------------------------------------------------------
    if time_col:
        for compound, sub_df in df.groupby("compound"):
            plt.figure(figsize=(6, 4))
            sns.lineplot(x=sub_df[time_col], y=sub_df["inhibition_%"], lw=2)
            plt.xlabel("Time (min)")
            plt.ylabel("% Inhibition")
            plt.title(f"{compound} - {target if target_col else ''}")
            plt.tight_layout()
            plt.savefig(output_dir / f"{compound}_inhibition_curve.png", dpi=150)
            plt.close()

    return {"df": result_df}

@register_task("pkpd_ranking", category="PBPK-PD",
               description="Rank compounds by predicted efficacy (e.g., Emax, AUEC, T>50%) and safety margins.")
def pkpd_ranking(config, df=None, **kwargs):
    """
    Rank compounds based on PK-PD simulation results.
    Optional weighting of multiple metrics (Emax, AUEC, T>50%, inhibition at Cmax, etc.).
    """
    if df is None:
        raise ValueError("Missing PK-PD analysis dataframe for ranking.")

    ranking_cfg = config.get("pkpd_ranking", {})
    output_dir = Path(ranking_cfg.get("output_dir", "outputs/pbpk/pkpd_ranking"))
    output_dir.mkdir(parents=True, exist_ok=True)

    # Choose which metric to rank by
    metric = ranking_cfg.get("metric", "Emax_mean")  # default to Emax_mean
    ascending = ranking_cfg.get("ascending", False)  # higher is better by default

    # Ensure the metric exists in df
    if metric not in df.columns:
        raise ValueError(f"Metric '{metric}' not found in PKPD dataframe. Available: {list(df.columns)}")

    # Optional weighting: combine multiple metrics
    weights = ranking_cfg.get("weights", None)  # dict of metric -> weight
    if weights:
        # Compute weighted score
        score = np.zeros(len(df))
        for m, w in weights.items():
            if m not in df.columns:
                raise ValueError(f"Metric '{m}' not in dataframe for weighting.")
            score += w * df[m].fillna(0).values
        df["ranking_score"] = score
        sort_col = "ranking_score"
    else:
        sort_col = metric

    # Sort compounds
    ranked_df = df.sort_values(by=sort_col, ascending=ascending).reset_index(drop=True)
    ranked_df["rank"] = ranked_df.index + 1

    # Save results
    ranked_path = output_dir / "pkpd_ranking.csv"
    ranked_df.to_csv(ranked_path, index=False)
    logger.info(f"PK–PD ranking saved to {ranked_path}")

    return {"df": ranked_df}

