import importlib
# import json

import pandas as pd
import numpy as np
# from joblib import Parallel, delayed
# from tqdm import tqdm

# from rdkit import Chem
# from rdkit.Chem import AllChem
# from rdkit.Chem import Descriptors
# from rdkit.Chem import QED

# from pathlib import Path

import torch
# from torch_geometric.loader import DataLoader
# from sklearn.model_selection import train_test_split
# from sklearn.preprocessing import StandardScaler, normalize
# from sklearn.decomposition import PCA

from pipeline.task_registry import register_task

from modules.utils.transforms import (ic50_to_pic50, ppb_to_logfu, vd_to_log, bioavailability_to_logit)

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def apply_saved_label_transform(series, metadata):
    transform = metadata.get("transform", "identity")

    values = pd.to_numeric(series, errors="coerce").astype(float)

    if transform == "ic50_to_pic50":
        values = values.apply(
            lambda x: ic50_to_pic50(x)
            if pd.notna(x) and x > 0 else np.nan
        )

    elif transform == "identity":
        clip = metadata.get("outlier_clip")
        if clip is not None:
            lower, upper = clip
            values.loc[(values < lower) | (values > upper)] = np.nan

    elif transform == "log":
        values.loc[values <= 0] = np.nan

        upper = metadata.get("upper_clip")
        if upper is not None:
            values.loc[values > upper] = upper

        mask = values.notna()
        values.loc[mask] = np.log(
            np.clip(values.loc[mask], 1e-6, None)
        )

    elif transform == "ppb_to_logfu":
        invalid = (
            values.notna()
            & (
                (values <= 0)
                | (values >= 100)
            )
        )
        values.loc[invalid] = np.nan

        mask = values.notna()
        values.loc[mask] = values.loc[mask].apply(ppb_to_logfu)

    elif transform == "log_vd":
        values.loc[values <= 0] = np.nan

        upper = metadata.get("upper_clip")
        if upper is not None:
            values.loc[values > upper] = upper

        mask = values.notna()
        values.loc[mask] = values.loc[mask].apply(vd_to_log)

    elif transform == "logit_f":
        invalid = (
            values.notna()
            & (
                (values <= 0)
                | (values >= 100)
            )
        )
        values.loc[invalid] = np.nan

        mask = values.notna()
        values.loc[mask] = values.loc[mask].apply(
            bioavailability_to_logit
        )

    elif transform == "log1p":
        invalid_mask = values.notna() & (~np.isfinite(values))
        values.loc[invalid_mask] = np.nan

        too_small = values.notna() & (values <= -0.999999)
        values.loc[too_small] = np.nan

        mask = values.notna()
        values.loc[mask] = np.log1p(values.loc[mask])

    else:
        raise ValueError(f"Unknown saved label transform: {transform}")

    return values


@register_task(
    "transform_labels_finetune",
    category="ADME",
    description="Apply pretrained ADME label transforms/scalers for fine-tuning."
)
def transform_labels_finetune(config, context):
    df = context["dataframe"].copy()

    checkpoint = context.get("adme_checkpoint")

    if checkpoint is None:
        raise RuntimeError(
            "transform_labels_finetune requires load_adme_checkpoint "
            "to run first."
        )

    pretrained_task_names = checkpoint["task_names"]
    requested_label_cols = config.get("label_cols", pretrained_task_names)

    if requested_label_cols != pretrained_task_names:
        raise ValueError(
            "Fine-tuning label_cols must exactly match checkpoint task_names "
            "and order.\n"
            f"YAML label_cols: {requested_label_cols}\n"
            f"Checkpoint task_names: {pretrained_task_names}"
        )

    scalers = checkpoint["label_scalers"]
    transform_metadata = checkpoint["label_transform_metadata"]

    for col in pretrained_task_names:
        if col not in df.columns:
            logger.warning(
                f"Fine-tuning CSV is missing column '{col}'. "
                "Adding it as NaN."
            )
            df[col] = np.nan

        transformed = apply_saved_label_transform(
            df[col],
            transform_metadata.get(col, {"transform": "identity"})
        )

        scaler = scalers.get(col)

        if scaler is None:
            logger.warning(
                f"No pretrained label scaler found for task '{col}'. "
                "Leaving transformed values unscaled."
            )
            df[col] = transformed
            continue

        mask = transformed.notna()

        if mask.sum() == 0:
            df[col] = np.nan
            continue

        scaled = scaler.transform(
            transformed.loc[mask].values.reshape(-1, 1)
        ).flatten()

        df[col] = np.nan
        df.loc[mask, col] = scaled

    context["dataframe"] = df
    context["label_scalers"] = scalers
    context["label_transform_metadata"] = transform_metadata
    context["task_names"] = pretrained_task_names

    label_counts = {
        col: int(df[col].notna().sum())
        for col in pretrained_task_names
    }

    logger.info(f"Fine-tuning label counts: {label_counts}")

    return context

@register_task(
    "load_adme_checkpoint",
    category="ADME",
    description="Load pretrained ADME model checkpoint and preprocessing artifacts."
)
def load_adme_checkpoint(config, context):
    checkpoint_path = config["checkpoint_path"]

    checkpoint = torch.load(
        checkpoint_path,
        map_location="cpu",
        weights_only=False,
    )

    required_keys = [
        "model_state_dict",
        "task_names",
        "task_groups",
        "params",
        "label_scalers",
        "label_transform_metadata",
        "global_feature_scaler",
    ]

    missing = [k for k in required_keys if k not in checkpoint]

    if missing:
        raise KeyError(
            f"Checkpoint is missing required keys for fine-tuning: {missing}"
        )

    context["adme_checkpoint"] = checkpoint
    context["pretrained_checkpoint_path"] = checkpoint_path
    context["label_scalers"] = checkpoint["label_scalers"]
    context["label_transform_metadata"] = checkpoint["label_transform_metadata"]
    context["task_names"] = checkpoint["task_names"]
    context["task_groups"] = checkpoint["task_groups"]

    logger.info(f"Loaded ADME checkpoint from {checkpoint_path}")
    logger.info(f"Checkpoint tasks: {checkpoint['task_names']}")

    return context

@register_task(
    "fine_tune_adme_model",
    category="ADME",
    description="Fine-tune pretrained multitask ADME model."
)
def fine_tune_adme_model(config, context):
    trainer_cfg = config["trainer"]

    module = importlib.import_module(trainer_cfg["module"])
    fine_tune_fn = getattr(module, trainer_cfg["function"])

    context["device"] = torch.device(
        config.get(
            "device",
            "cuda" if torch.cuda.is_available() else "cpu"
        )
    )

    result = fine_tune_fn(
        context=context,
        config=config
    )

    trained_model = result.get("model")

    context.update({
        "trained_model": trained_model,
        "training_summary": result,
    })

    return context