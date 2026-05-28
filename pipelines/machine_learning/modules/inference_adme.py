import importlib
import joblib

from pathlib import Path
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import QED

import torch
from torch_geometric.loader import DataLoader
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, normalize
from sklearn.decomposition import PCA

from pipeline.task_registry import register_task

from modules.utils.plotting import plot_label_histograms
from modules.utils..conversions import logfu_to_ppb, log_to_vd, logit_to_bioavailability
from modules.utils.splits import scaffold_split

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_task("load_smiles_inference_dataset", category="ADME")
def load_smiles_inference_dataset(config, context):
    df = pd.read_csv(config["csv_path"])

    smiles_col = config.get("smiles_col", "smiles")

    if smiles_col not in df.columns:
        raise ValueError(f"{smiles_col} not found in CSV")

    context["dataframe"] = df
    context["smiles_col"] = smiles_col

    return context

@register_task("featurise_smiles_inference", category="ADME")
def featurise_smiles_inference(config, context):
    df = context["dataframe"]
    smiles_col = context["smiles_col"]

    feat_cfg = config["featuriser"]
    module = importlib.import_module(feat_cfg["module"])
    featuriser = getattr(module, feat_cfg["function"])

    if hasattr(module, "prepare_features"):
        module.prepare_features(df, smiles_col=smiles_col)

    graphs = []
    valid_indices = []

    for i, smi in enumerate(tqdm(df[smiles_col], desc="Featurising")):

        # -------------------------
        # Skip missing SMILES
        # -------------------------
        if pd.isna(smi):
            logger.warning(f"Skipping row {i}: SMILES is NaN")
            continue

        smi = str(smi).strip()

        if smi == "":
            logger.warning(f"Skipping row {i}: SMILES is empty")
            continue

        # -------------------------
        # Featurise safely
        # -------------------------
        try:
            g = featuriser(smi, label=None, idx=i)

            if g is not None:
                graphs.append(g)
                valid_indices.append(i)
            else:
                logger.warning(f"Skipping row {i}: featuriser returned None")

        except Exception as e:
            logger.warning(f"Failed to featurise row {i} ({smi}): {e}")
            continue

    if not graphs:
        raise ValueError("No valid molecules for inference.")

    first_graph = graphs[0]

    context.update({
        "graphs": graphs,
        "valid_indices": valid_indices,
        "input_dim": first_graph.x.shape[1],
        "edge_dim": first_graph.edge_attr.shape[1],
        "global_feat_dim": first_graph.global_features.shape[-1],
        "fp_dim": first_graph.fp.shape[-1],
    })

    return context

@register_task("load_trained_adme_model", category="ADME")
def load_trained_adme_model(config, context):
    model_path = config["model_path"]

    checkpoint = torch.load(model_path, map_location="cpu")

    # get scalers 

    context["label_scalers"] = checkpoint.get("label_scalers")
    context["label_transform_metadata"] = checkpoint.get("label_transform_metadata")

    # -------------------------
    # Extract metadata
    # -------------------------
    input_dim = checkpoint["input_dim"]
    edge_dim = checkpoint["edge_dim"]
    global_feat_dim = checkpoint["global_feat_dim"]
    fp_dim = checkpoint["fp_dim"]
    num_tasks = checkpoint["num_tasks"]
    params = checkpoint.get("params", {})
    task_names = checkpoint.get("task_names", [])

    context["task_names"] = task_names

    # -------------------------
    # Task groups (REQUIRED)
    # -------------------------
    task_groups_cfg = config.get("task_groups")

    if not task_groups_cfg:
        raise ValueError("task_groups must be provided in inference config")

    # Convert names → indices
    task_groups = {}
    for group, names in task_groups_cfg.items():
        indices = [task_names.index(n) for n in names]
        task_groups[group] = indices

    context["task_groups"] = task_groups

    # -------------------------
    # Build model
    # -------------------------
    ModelClass = context["model_class"]

    model = ModelClass(
        input_dim=input_dim,
        edge_dim=edge_dim,
        global_feat_dim=global_feat_dim,
        fp_dim=fp_dim,
        num_tasks=num_tasks,
        task_groups=task_groups,
        **{k: v for k, v in params.items() if k != "lr"}
    )

    # -------------------------
    # Load weights
    # -------------------------
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()

    context["model"] = model

    return {"model": model}

@register_task("run_adme_inference", category="ADME")
def run_adme_inference(config, context):
    graphs = context["graphs"]
    model = context["model"]
    scalers = context.get("label_scalers")
    task_names = context["task_names"]

    batch_size = config.get("batch_size", 64)
    loader = DataLoader(graphs, batch_size=batch_size, shuffle=False)

    preds = []

    # -------------------------
    # RUN MODEL
    # -------------------------
    with torch.no_grad():
        for batch in tqdm(loader, desc="Inference"):
            batch = batch.to("cpu")
            out = model(batch)
            preds.append(out.cpu().numpy())

    preds = np.vstack(preds)
    
    # -------------------------
    # INVERSE SCALING + TRANSFORMS
    # -------------------------
    transform_metadata = context.get("label_transform_metadata")

    final_preds = np.zeros_like(preds)

    for i, task in enumerate(task_names):

        values = preds[:, i]

        # ---- Step 1: inverse scaling ----
        if scalers is not None and task in scalers:
            scaler = scalers[task]
            values = scaler.inverse_transform(values.reshape(-1, 1)).flatten()

        # ---- Step 2: inverse transform ----
        if transform_metadata and task in transform_metadata:

            meta = transform_metadata[task]
            transform_type = meta.get("transform")

            if transform_type == "log":
                values = np.exp(values)

            elif transform_type == "log1p":
                values = np.expm1(values)

            elif transform_type == "ic50_to_pic50":
                # pIC50 -> IC50 nM
                values = 10 ** (-values) * 1e9

            elif transform_type == "ppb_to_logfu":
                values = logfu_to_ppb(values)

            elif transform_type == "log_vd":
                values = log_to_vd(values)

            elif transform_type == "logit_f":
                values = logit_to_bioavailability(values)

            elif transform_type == "identity":
                pass

            else:
                logger.warning(f"Unknown transform {transform_type} for task {task}")

        final_preds[:, i] = values

    context["predictions"] = final_preds
    context["task_names"] = task_names

    return context

@register_task("save_adme_predictions", category="ADME")
def save_adme_predictions(config, context):
    df = context["dataframe"].copy()
    preds = context["predictions"]
    task_names = context["task_names"]
    valid_indices = context["valid_indices"]

    pred_df = pd.DataFrame(preds, columns=task_names)

    # align with original dataframe
    for col in task_names:
        df[col + "_pred"] = np.nan

    df.loc[
        valid_indices,
        [col + "_pred" for col in task_names]
    ] = pred_df.values

    output_path = Path(config["output_path"])

    # Create directory if needed
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(output_path, index=False)

    return {"output_path": str(output_path)}