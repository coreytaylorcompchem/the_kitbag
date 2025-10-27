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

def get_available_pbpk_models():
    """
    Returns a mapping of available PBPK model names to their corresponding
    file paths, classes, and featurisers.
    """
    # Base models directory relative to pipeline root
    model_dir = Path(__file__).resolve().parent.parent / "models"

    models = {
        "clint_model": {
            "path": model_dir / "clint" / "hep_clint_gin.pt",
            "class": ClintModel,
            "featuriser": mol_to_graph_clint,
        },
        "logbb_model": {
            "path": model_dir / "logbb" / "log_bb_gin.pt",
            "class": LogBBModel,
            "featuriser": mol_to_graph_logbb,
        },
        "fu_p_model": {
            "path": model_dir / "ppb_fu" / "ppb_f_u_gin.pt",
            "class": PPBFuModel,
            "featuriser": mol_to_graph_ppbfu,
        },
    }

    return models

# ----------------------------------------------------------
# 1. Parameter Prediction (minimal PoC version)
# ----------------------------------------------------------
@register_task("pbpk_parameter_prediction", category="PBPK", description="Predict PBPK parameters from SMILES")
def pbpk_parameter_prediction(config, df=None):
    """Predict PBPK parameters from SMILES."""
    pbpk_cfg = config.get("pbpk_parameter_prediction", {})
    predictors_cfg = pbpk_cfg.get("predictors", {})
    output_dir = Path(pbpk_cfg.get("output_dir", "outputs/pbpk/params"))
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"[Config check] Predictors in YAML: {predictors_cfg}")

    all_models = get_available_pbpk_models()
    selected_models = {
        name: all_models[name]
        for name, enabled in predictors_cfg.items()
        if isinstance(enabled, bool) and enabled and name in all_models
    }

    if not selected_models:
        raise ValueError(
            "No PBPK models enabled in YAML under 'pbpk_parameter_prediction.predictors'. "
            f"Current config: {predictors_cfg}"
        )

    logger.info(f"Enabled PBPK models: {list(selected_models.keys())}")

    input_file = config.get("input_file")
    if not input_file:
        raise ValueError("Missing 'input_file' in config.")

    df = pd.read_csv(input_file)
    if "smiles" not in df.columns:
        raise ValueError("Input file must contain a 'smiles' column.")

    # DEBUGGING: any weird / failed SMILES
    for idx, smi in enumerate(df["smiles"]):
        s = str(smi)
        s_stripped = s.strip()
        logger.debug(f"Row {idx}: original='{s}' (len={len(s)}) | stripped='{s_stripped}' (len={len(s_stripped)}) | repr={repr(s)}")

    # Simple filter for poor SMILES
    df["smiles"] = df["smiles"].astype(str).str.strip()
    df = df[df["smiles"].notna() & (df["smiles"] != "")].copy()

    result_df = df[["name", "smiles"]].copy()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    logger.info(f"Using device: {device}")

    # Run predictions for each enabled model
    for model_name, info in selected_models.items():
        model_path = info["path"]
        model_class = info["class"]
        featuriser = info["featuriser"]

        if not model_path.exists():
            raise FileNotFoundError(f"Model file not found at: {model_path}")

        logger.info(f"Loading model '{model_name}' from {model_path}")
        checkpoint = torch.load(model_path, map_location=device)

        sig = inspect.signature(model_class.__init__)
        init_params = sig.parameters

        model_args = {}
        for arg in list(init_params.keys())[1:]:
            if arg in checkpoint:
                model_args[arg] = checkpoint[arg]
            elif arg in checkpoint.get("model_hparams", {}):
                model_args[arg] = checkpoint["model_hparams"][arg]
            elif arg == "input_dim":
                if "input_dim" in checkpoint:
                    model_args[arg] = checkpoint["input_dim"]
                elif "input_dim" in checkpoint.get("model_hparams", {}):
                    model_args[arg] = checkpoint["model_hparams"]["input_dim"]
                else:
                    example_graph = featuriser("CCO")
                    if example_graph is not None and hasattr(example_graph, "x"):
                        model_args[arg] = example_graph.x.shape[1]
                    else:
                        model_args[arg] = 0
            elif arg == "edge_dim":
                if "edge_dim" in checkpoint:
                    model_args[arg] = checkpoint["edge_dim"]
                elif "edge_dim" in checkpoint.get("model_hparams", {}):
                    model_args[arg] = checkpoint["model_hparams"]["edge_dim"]
                else:
                    example_graph = featuriser("CCO")
                    if example_graph is not None and hasattr(example_graph, "edge_attr"):
                        model_args[arg] = example_graph.edge_attr.shape[1]
                    else:
                        model_args[arg] = 0
            elif arg == "global_feat_dim":
                if "global_feat_dim" in checkpoint:
                    model_args[arg] = checkpoint["global_feat_dim"]
                elif "global_feat_dim" in checkpoint.get("model_hparams", {}):
                    model_args[arg] = checkpoint["model_hparams"]["global_feat_dim"]
                else:
                    example_graph = featuriser("CCO")
                    if example_graph is not None and hasattr(example_graph, "global_features"):
                        model_args[arg] = example_graph.global_features.shape[0]
                    else:
                        model_args[arg] = 0
            elif arg == "hidden_dim":
                model_args[arg] = checkpoint.get("hidden_dim", 128)
            else:
                # fallback for optional args with default values
                if init_params[arg].default is not inspect._empty:
                    model_args[arg] = init_params[arg].default
                else:
                    logger.warning(f"Could not infer value for '{arg}' when loading {model_name}")

        logger.info(f"Instantiating {model_name} with args: {model_args}")
        model = model_class(**model_args)
        model.load_state_dict(checkpoint["model_state_dict"])
        model.to(device)
        model.eval()

        #  Run inference for DL models
        preds = []
        failed_smiles = []

        for smi in df["smiles"]:
            graph_data = featuriser(smi)
            if graph_data is None:
                preds.append(None)
                failed_smiles.append(smi)
                continue

            loader = DataLoader([graph_data], batch_size=1, shuffle=False)
            with torch.no_grad():
                for batch in loader:
                    batch = batch.to(device)
                    output = model(batch)
                    preds.append(output.item())

        if failed_smiles:
            logger.info(f"Featurisation failed for {len(failed_smiles)} SMILES: {failed_smiles}")

        result_df[model_name.replace("_model", "")] = preds

    out_path = output_dir / "pbpk_params.csv"
    result_df.to_csv(out_path, index=False)
    logger.info(f"Parameter predictions saved to {out_path}")

    return {"df": result_df}

# ----------------------------------------------------------
# 2. PBPK model assembly
# ----------------------------------------------------------

@register_task("pbpk_model_assembly", category="PBPK", description="Assemble PBPK model")
def pbpk_model_assembly(config, df=None):
    """
    Assemble the PBPK model by combining physiology constants and ML-predicted parameters.
    Returns both physiology and the parameter dataframe for downstream simulation.
    """

    assembly_cfg = config.get("pbpk_model_assembly", {})
    physiology_file = assembly_cfg.get("physiology_file")
    output_dir = Path(assembly_cfg.get("output_dir", "outputs/pbpk/models"))
    output_dir.mkdir(parents=True, exist_ok=True)

    if not physiology_file:
        raise ValueError("Missing 'physiology_file' in pbpk_model_assembly configuration.")

    physiology_path = Path(os.path.expandvars(os.path.expanduser(physiology_file))).resolve()
    if not physiology_path.exists():
        raise FileNotFoundError(f"Physiology file not found at: {physiology_path}")

    logger.info(f"Using physiology file: {physiology_path}")
    with open(physiology_path, "r") as f:
        physiology = yaml.safe_load(f)

    brain_phys = physiology.get("physiology", {}).get("brain", {})
    logger.debug(f"Loaded brain physiology keys: {list(brain_phys.keys())}")

    # Save combined physiology + model parameters
    assembled_path = output_dir / "pbpk_assembled_model.json"
    with open(assembled_path, "w") as f:
        json.dump(physiology, f, indent=2)

    if df is not None:
        params_path = output_dir / "pbpk_params_assembled.csv"
        df.to_csv(params_path, index=False)
        logger.info(f"Saved predicted PBPK parameters to {params_path}")
    else:
        logger.warning("No DataFrame (df) provided; downstream simulation may fail.")

    logger.info(f"PBPK model assembled and saved to {output_dir}")

    result = {
        "physiology": physiology.get("physiology", physiology),
        "df": df
    }

    logger.debug(f"Returning assembly context keys: {list(result.keys())}")
    return result


# ----------------------------------------------------------
# 3. PBPK Simulation
# ----------------------------------------------------------

# Helper: ensure fu_p is a proportion
def fix_fu_p(fu_p_percentage):
    return min(max(fu_p_percentage / 100, 0), 1)

def generate_virtual_subjects(base_phys, base_params, n_subjects, variability):
    """
    Generates virtual subjects by applying log-normal variability to PBPK parameters.
    Returns a list of dicts with subject-specific parameters.
    """
    subjects = []
    flows = base_phys.get('flows', {})
    for _ in range(n_subjects):
        subj = {}
        # Plasma & brain volumes
        sigma_Vp = max(base_phys['plasma'].get('V_p', 3.0) * variability.get('Vp', 0), 0)
        subj['V_p'] = np.random.normal(base_phys['plasma'].get('V_p', 3.0), sigma_Vp)

        sigma_Vb = max(base_phys['brain'].get('V_b', 1.0) * variability.get('Vb', 0), 0)
        subj['V_b'] = np.random.normal(base_phys['brain'].get('V_b', 1.0), sigma_Vb)

        # Brain flow & PS
        sigma_Q = max(flows.get('Q_brain', 0.75) * variability.get('Q_brain', 0), 0)
        subj['Q_brain'] = np.random.normal(flows.get('Q_brain', 0.75), sigma_Q)

        sigma_PS = max(base_phys['brain'].get('PS_BBB', 0.05) * variability.get('PS_BBB', 0), 0)
        subj['PS_BBB'] = np.random.normal(base_phys['brain'].get('PS_BBB', 0.05), sigma_PS)

        # F_u plasma
        sigma_fu = max(base_params.get('fu_p', 0.1) * variability.get('fu_p', 0), 0)
        subj['fu_p'] = np.random.normal(base_params.get('fu_p', 0.1), sigma_fu)

        subjects.append(subj)
    return subjects

@register_task("pbpk_simulation", category="PBPK - brain:plasma",
               description="Run 2-compartment plasma–brain PBPK simulation with selected exchange model.")
def pbpk_simulation(config, df=None, physiology=None, **kwargs):
    """
    Runs a 2-compartment PBPK simulation (plasma-brain) using provided parameters.
    Combines ML-predicted PK parameters with physiology constants and generates virtual subjects.
    Aggregates results across subjects and plots mean ± SD curves and 95% CI.
    """
    if df is None:
        raise ValueError("Missing 'df' (predicted PBPK parameters) for PBPK simulation.")
    if physiology is None:
        raise ValueError("Missing 'physiology' data for PBPK simulation.")

    sim_cfg = config.get("pbpk_simulation", {})
    t_max = sim_cfg.get("simulation_time", 240)
    dose = sim_cfg.get("dose_mg", 10.0)
    exchange_model = sim_cfg.get("exchange_model", "hybrid").lower()
    n_subjects = sim_cfg.get("n_subjects", 5)
    variability = sim_cfg.get("variability", {"Vp":0.1, "Vb":0.1, "Q_brain":0.1, "PS_BBB":0.1, "fu_p":0.1})

    outdir = Path(sim_cfg.get("output_dir", "outputs/pbpk/simulations"))
    outdir.mkdir(parents=True, exist_ok=True)

    indiv_dir = outdir / "individual_plots"
    indiv_dir.mkdir(exist_ok=True)

    phys_root = physiology.get("physiology", physiology)
    
    results = []
    aggregate_results = []

    z95 = 1.96  # Z-score CIs; 95% CI multiplier

    # Per-compound loop

    for _, row in tqdm(df.iterrows(), total=len(df), desc="Compounds", unit="compound"):
        subjects = generate_virtual_subjects(phys_root, row.to_dict(), n_subjects, variability)
        compound_name = row.get("compound", row.get("name", "unknown"))

        Cp_traces = []
        Cb_traces = []

        # Per-subject loop for population calculations

        for subj_idx, subj in tqdm(enumerate(subjects), total=n_subjects,
                               desc=f"{compound_name} subjects", leave=False, unit="subject"):
            Cp0, Cb0 = 0.0, 0.0
            fu_p = max(subj.get("fu_p", 0.1), 1e-6)

            Kp = 10 ** row["logbb"]
            CLint = row["clint"] / 60.0
            PS = subj.get("PS_BBB", phys_root.get("brain", {}).get("PS_BBB", 0.05))
            k_efflux = subj.get("k_efflux", phys_root.get("brain", {}).get("k_efflux", 0.001))
            V_p = subj.get("V_p", phys_root.get("plasma", {}).get("V_p", 3.0))
            V_b = subj.get("V_b", phys_root.get("brain", {}).get("V_b", 1.0))
            Q_brain = subj.get("Q_brain", phys_root.get("flows", {}).get("Q_brain", 0.75))

            if exchange_model == "flow_limited":
                def brain_flux(Cp, Cb):
                    return Q_brain * (Cp - Cb / Kp)
            elif exchange_model == "permeability_limited":
                def brain_flux(Cp, Cb):
                    return PS * (Cp - Cb / Kp)
            else:
                def brain_flux(Cp, Cb):
                    eff_flux = (Q_brain * PS) / (Q_brain + PS)
                    return eff_flux * (Cp - Cb / Kp)

            def ode(t, y):
                Cp, Cb = y
                input_term = (dose / V_p) if t < 1 else 0
                J = brain_flux(Cp, Cb)
                dCp = input_term - (J / V_p) - fu_p * CLint * Cp
                dCb = (J / V_b) - k_efflux * Cb
                return [dCp, dCb]

            t_eval = np.linspace(0, t_max, 200)
            sol = solve_ivp(ode, [0, t_max], [Cp0, Cb0], t_eval=t_eval)
            Cp, Cb = sol.y

            Cp_traces.append(Cp)
            Cb_traces.append(Cb)

            auc_p = np.trapz(Cp, sol.t)
            auc_b = np.trapz(Cb, sol.t)
            cmax_p = Cp.max()
            cmax_b = Cb.max()

            results.append({
                "compound": compound_name,
                "subject": subj_idx,
                "exchange_model": exchange_model,
                "AUCp": auc_p,
                "AUCb": auc_b,
                "Cmax_plasma": cmax_p,
                "Cmax_brain": cmax_b,
                "Kp_brain": auc_b / auc_p if auc_p > 0 else np.nan,
                "V_p": V_p,
                "V_b": V_b,
                "Q_brain": Q_brain,
                "PS_BBB": PS,
                "k_efflux": k_efflux,
            })

            # Set palette for all plots
            palette = sns.color_palette("Dark2")  # can add any seaborn palette or list of colors
            plt.rcParams['axes.prop_cycle'] = cycler(color=palette)

            plt.figure(figsize=(5, 3))
            plt.plot(t_eval, Cp, label="Plasma", lw=1.8)
            plt.plot(t_eval, Cb, label="Brain", lw=1.8)
            plt.legend()
            plt.xlabel("Time (min)")
            plt.ylabel("Concentration (a.u.)")
            plt.title(f"{compound_name} - Subject {subj_idx} ({exchange_model})")
            plt.tight_layout()
            plt.savefig(indiv_dir / f"{compound_name}_subject{subj_idx}_{exchange_model}.png", dpi=150)
            plt.close()

        # Aggregate all metrics
        comp_df = pd.DataFrame(results[-n_subjects:])

        def ci95(x):
            return z95 * sem(x) if len(x) > 1 else np.nan

        aggregate_results.append({
            "compound": compound_name,
            "exchange_model": exchange_model,
            "AUCp_mean": comp_df["AUCp"].mean(),
            "AUCp_std": comp_df["AUCp"].std(),
            "AUCp_95CI": ci95(comp_df["AUCp"]),
            "AUCb_mean": comp_df["AUCb"].mean(),
            "AUCb_std": comp_df["AUCb"].std(),
            "AUCb_95CI": ci95(comp_df["AUCb"]),
            "Cmax_plasma_mean": comp_df["Cmax_plasma"].mean(),
            "Cmax_plasma_std": comp_df["Cmax_plasma"].std(),
            "Cmax_plasma_95CI": ci95(comp_df["Cmax_plasma"]),
            "Cmax_brain_mean": comp_df["Cmax_brain"].mean(),
            "Cmax_brain_std": comp_df["Cmax_brain"].std(),
            "Cmax_brain_95CI": ci95(comp_df["Cmax_brain"]),
            "Kp_brain_mean": comp_df["Kp_brain"].mean(),
            "Kp_brain_std": comp_df["Kp_brain"].std(),
            "Kp_brain_95CI": ci95(comp_df["Kp_brain"]),
        })

        # Aggregate plot (mean ± SD)
        Cp_traces = np.array(Cp_traces)
        Cb_traces = np.array(Cb_traces)

        Cp_mean = Cp_traces.mean(axis=0)
        Cp_std = Cp_traces.std(axis=0)
        Cb_mean = Cb_traces.mean(axis=0)
        Cb_std = Cb_traces.std(axis=0)

        Cp_95CI = ci95(Cp_traces)
        Cb_95CI = ci95(Cb_traces)

        # Define traces and labels
        traces = [
            (Cp_mean, Cp_std, Cp_95CI, "Plasma"),
            (Cb_mean, Cb_std, Cb_95CI, "Brain")
        ]

        fig, ax = plt.subplots(figsize=(6, 4))
        
        prop_cycle = plt.rcParams['axes.prop_cycle']  # get default property cycle
        colors = [d['color'] for d in prop_cycle]

        # Example usage for multiple traces
        for i, (y_mean, y_std, y_ci, label) in enumerate(traces):
            color = colors[i % len(colors)]
            ax.plot(t_eval, y_mean, lw=2, label=f"{label} Mean", color=color)
            ax.fill_between(t_eval, y_mean - y_std, y_mean + y_std, color=color, alpha=0.3, label=f"{label} ±1 SD")
            ax.fill_between(t_eval, y_mean - y_ci, y_mean + y_ci, color=color, alpha=0.2, label=f"{label} 95% CI")
    
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("Concentration (a.u.)")
        ax.set_title(f"{compound_name} Aggregate ({exchange_model})")
        ax.legend()
        plt.tight_layout()
        plt.savefig(outdir / f"{compound_name}_aggregate_{exchange_model}.png", dpi=150)
        plt.close()

    # Save results to csv
    result_df = pd.DataFrame(results)
    agg_df = pd.DataFrame(aggregate_results)
    result_df.to_csv(outdir / "pbpk_simulation_results.csv", index=False)
    agg_df.to_csv(outdir / "pbpk_simulation_aggregate.csv", index=False)

    logger.info(f"Per-subject results saved to: {outdir / 'pbpk_simulation_results.csv'}")
    logger.info(f"Aggregate results saved to: {outdir / 'pbpk_simulation_aggregate.csv'}")
    logger.info(f"Individual plots saved to: {indiv_dir}")

    return {
        "df": result_df,
        "df_aggregate": agg_df,
        "physiology": physiology
    }


# ----------------------------------------------------------
# 4. Analysis / Derived Metrics
# ----------------------------------------------------------
@register_task("pbpk_analysis", category="PBPK", description="Analyse PBPK simulation results.")
def pbpk_analysis(config, df=None, physiology=None, **kwargs):
    """
    Perform post-simulation analysis of PBPK results.
    Computes metrics: auc_plasma, auc_brain, kp_brain, kp_uu, cmax, cmax_ratio.
    Accepts simulation output DataFrame (df) and optional physiology data.
    """

    if df is None:
        raise ValueError("Missing 'df' containing PBPK simulation results for analysis.")

    analysis_cfg = config.get("pbpk_analysis", {})
    outdir = Path(analysis_cfg.get("output_dir", "outputs/pbpk/analysis"))
    outdir.mkdir(parents=True, exist_ok=True)

    # Compute metrics
    metrics_df = pd.DataFrame()
    metrics_df["compound"] = df["compound"]
    metrics_df["auc_plasma"] = df["AUCp"]
    metrics_df["auc_brain"] = df["AUCb"]
    metrics_df["kp_brain"] = df["Kp_brain"]

    # Compute unbound brain:plasma ratio if fu_p is present
    if "fu_p" in df.columns:
        metrics_df["kp_uu"] = metrics_df["kp_brain"] / df["fu_p"].apply(lambda x: max(x/100.0, 1e-6))
    else:
        metrics_df["kp_uu"] = metrics_df["kp_brain"]
    # Cmax and Cmax ratio
    metrics_df["cmax_plasma"] = df["Cmax_plasma"]
    metrics_df["cmax_brain"] = df["Cmax_brain"]
    metrics_df["cmax_ratio"] = df["Cmax_brain"] / df["Cmax_plasma"].replace(0, np.nan)

    summary_path = outdir / "pbpk_summary.csv"
    metrics_df.to_csv(summary_path, index=False)
    logger.info(f"Saved PBPK metrics summary to {summary_path}")

    if physiology:
        physio_path = outdir / "physiology_used.json"
        with open(physio_path, "w") as f:
            json.dump(physiology, f, indent=2)
        logger.debug(f"Saved physiology snapshot to {physio_path}")

    return {"df": metrics_df, "physiology": physiology}

