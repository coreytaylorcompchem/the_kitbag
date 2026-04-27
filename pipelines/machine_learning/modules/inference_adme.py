import importlib
import joblib

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

    # ✅ THIS IS THE MISSING LINE
    if hasattr(module, "prepare_features"):
        module.prepare_features(df, smiles_col=smiles_col)

    graphs = []
    valid_indices = []

    for i, smi in enumerate(tqdm(df[smiles_col], desc="Featurising")):
        g = featuriser(smi, label=None, idx=i)

        if g is not None:
            graphs.append(g)
            valid_indices.append(i)

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
    # device = context["device"]
    model = context["model"]
    scalers = context.get("label_scalers")
    task_names = context["task_names"]

    if scalers is not None:
        preds_rescaled = np.zeros_like(preds)

        for i, task in enumerate(task_names):
            scaler = scalers[task]
            preds_rescaled[:, i] = scaler.inverse_transform(
                preds[:, i].reshape(-1, 1)
            ).flatten()

        final_preds = preds_rescaled
    else:
        logger.warning("Skipping inverse scaling — returning normalised predictions.")
        final_preds = preds

    batch_size = config.get("batch_size", 64)

    loader = DataLoader(graphs, batch_size=batch_size, shuffle=False)

    preds = []

    with torch.no_grad():
        for batch in tqdm(loader, desc="Inference"):
            batch = batch.to('cpu')
            out = model(batch)
            preds.append(out.cpu().numpy())

    preds = np.vstack(preds)

    # -------------------------
    # INVERSE SCALING
    # -------------------------
    task_names = list(scalers.keys())

    preds_rescaled = np.zeros_like(preds)

    for i, task in enumerate(task_names):
        scaler = scalers[task]
        preds_rescaled[:, i] = scaler.inverse_transform(
            preds[:, i].reshape(-1, 1)
        ).flatten()

    context["predictions"] = preds_rescaled
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

    df.loc[valid_indices, [col + "_pred" for col in task_names]] = pred_df.values

    output_path = config["output_path"]
    df.to_csv(output_path, index=False)

    return {"output_path": output_path}