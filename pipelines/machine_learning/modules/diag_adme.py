
import os
import math
import inspect
import joblib
from collections import OrderedDict

import textwrap
import matplotlib.pyplot as plt

import importlib
import json

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
from modules.utils.transforms import (ic50_to_pic50, ppb_to_logfu, vd_to_log, bioavailability_to_logit)

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def _ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path


def _safe_float(x):
    if x is None:
        return np.nan

    try:
        x = float(x)
    except Exception:
        return np.nan

    if not np.isfinite(x):
        return np.nan

    return x


def _get_state_dict_from_checkpoint(checkpoint):
    if isinstance(checkpoint, OrderedDict):
        return checkpoint

    if not isinstance(checkpoint, dict):
        raise TypeError(
            "Checkpoint must be a dict or OrderedDict containing model weights."
        )

    candidate_keys = [
        "model_state_dict",
        "state_dict",
        "model",
        "model_state",
        "net",
    ]

    for key in candidate_keys:
        value = checkpoint.get(key)
        if isinstance(value, (dict, OrderedDict)):
            return value

    tensor_like_values = [
        value for value in checkpoint.values()
        if torch.is_tensor(value)
    ]

    if tensor_like_values:
        return checkpoint

    raise KeyError(
        "Could not find model state dict in checkpoint. Expected one of: "
        f"{candidate_keys}, or a raw state_dict."
    )


def _strip_module_prefix_from_state_dict(state_dict):
    cleaned = OrderedDict()

    for key, value in state_dict.items():
        if key.startswith("module."):
            cleaned[key[len("module."):]] = value
        else:
            cleaned[key] = value

    return cleaned


def _checkpoint_get(checkpoint, keys, default=None):
    if not isinstance(checkpoint, dict):
        return default

    for key in keys:
        if key in checkpoint:
            return checkpoint[key]

    return default


@register_task(
    "load_adme_checkpoint",
    category="ADME",
    description="Load a trained ADME checkpoint plus stored preprocessing state."
)
def load_adme_checkpoint(config, context):
    checkpoint_path = config["checkpoint_path"]
    device = torch.device(config.get("device", "cpu"))

    logger.info(f"Loading ADME checkpoint from: {checkpoint_path}")

    checkpoint = torch.load(
        checkpoint_path,
        map_location=device,
        weights_only=False,
    )

    label_scalers = _checkpoint_get(
        checkpoint,
        keys=["label_scalers", "scalers", "target_scalers"],
        default=None,
    )

    transform_metadata = _checkpoint_get(
        checkpoint,
        keys=[
            "label_transform_metadata",
            "transform_metadata",
            "target_transform_metadata",
        ],
        default=None,
    )

    global_feature_scaler = _checkpoint_get(
        checkpoint,
        # keys=["global_feature_scaler", "feature_scaler"],
        keys=["feature_scaler"],
        default=None,
    )

    feature_state = _checkpoint_get(
        checkpoint,
        keys=["feature_state"],
        default=None,
    )

    scaler_path = config.get("scaler_path")
    transform_metadata_path = config.get("transform_metadata_path")

    if label_scalers is None and scaler_path:
        logger.info(f"Loading label scalers from: {scaler_path}")
        label_scalers = joblib.load(scaler_path)

    if transform_metadata is None and transform_metadata_path:
        logger.info(f"Loading transform metadata from: {transform_metadata_path}")
        with open(transform_metadata_path, "r") as f:
            transform_metadata = json.load(f)

    if label_scalers is None:
        raise RuntimeError(
            "No label scalers found in checkpoint and no scaler_path was provided."
        )

    if transform_metadata is None:
        raise RuntimeError(
            "No label transform metadata found in checkpoint and no "
            "transform_metadata_path was provided."
        )

    context["adme_checkpoint"] = checkpoint
    context["adme_checkpoint_path"] = checkpoint_path
    context["checkpoint_label_scalers"] = label_scalers
    context["checkpoint_label_transform_metadata"] = transform_metadata

    if global_feature_scaler is not None:
        context["global_feature_scaler"] = global_feature_scaler

    if feature_state is not None:
        context["checkpoint_feature_state"] = feature_state

    logger.info("Loaded ADME checkpoint preprocessing state")
    logger.info(f"Number of checkpoint label scalers: {len(label_scalers)}")
    logger.info(f"Number of checkpoint transform metadata entries: {len(transform_metadata)}")

    return context


def _apply_checkpoint_label_transform(series, metadata):
    series = pd.to_numeric(series, errors="coerce")

    transform_name = metadata.get("transform", "identity")

    if transform_name in {None, "none", "identity"}:
        out = series.copy()

        outlier_clip = metadata.get("outlier_clip")
        if outlier_clip is not None:
            lower, upper = outlier_clip
            mask = (out < lower) | (out > upper)
            out.loc[mask] = np.nan

        return out

    if transform_name in {"ic50_to_pic50", "ec50_to_pec50"}:
        return series.apply(
            lambda x: ic50_to_pic50(x)
            if pd.notna(x) and x > 0 else np.nan
        )

    if transform_name == "log":
        out = series.copy()
        out.loc[out <= 0] = np.nan

        upper = metadata.get("upper_clip")
        if upper is not None:
            out.loc[out > upper] = upper

        mask = out.notna()
        out.loc[mask] = np.log(np.clip(out.loc[mask], 1e-6, None))

        return out

    if transform_name == "log1p":
        out = series.copy()

        invalid_mask = out.notna() & (~np.isfinite(out))
        out.loc[invalid_mask] = np.nan

        too_small = out.notna() & (out <= -0.999999)
        out.loc[too_small] = np.nan

        mask = out.notna()
        out.loc[mask] = np.log1p(out.loc[mask])

        return out

    if transform_name == "ppb_to_logfu":
        out = series.copy()

        invalid = (
            out.notna() &
            (
                (out <= 0) |
                (out >= 100)
            )
        )

        out.loc[invalid] = np.nan

        mask = out.notna()
        out.loc[mask] = out.loc[mask].apply(ppb_to_logfu)

        return out

    if transform_name == "log_vd":
        out = series.copy()

        out.loc[out <= 0] = np.nan

        upper = metadata.get("upper_clip")
        if upper is not None:
            out.loc[out > upper] = upper

        mask = out.notna()
        out.loc[mask] = out.loc[mask].apply(vd_to_log)

        return out

    if transform_name == "logit_f":
        out = series.copy()

        invalid = (
            out.notna() &
            (
                (out <= 0) |
                (out >= 100)
            )
        )

        out.loc[invalid] = np.nan

        mask = out.notna()
        out.loc[mask] = out.loc[mask].apply(bioavailability_to_logit)

        return out

    raise ValueError(
        f"Unsupported checkpoint label transform: '{transform_name}'"
    )


@register_task(
    "transform_labels_from_checkpoint",
    category="ADME",
    description="Apply checkpoint label transforms and checkpoint scalers without refitting."
)
def transform_labels_from_checkpoint(config, context):
    df = context["dataframe"].copy()

    label_cols = config.get("label_cols")
    if label_cols is None:
        label_cols = context.get("task_names")

    if label_cols is None:
        raise ValueError(
            "transform_labels_from_checkpoint requires label_cols in YAML "
            "or context['task_names']."
        )

    label_scalers = context.get("checkpoint_label_scalers")
    transform_metadata = context.get("checkpoint_label_transform_metadata")

    if label_scalers is None:
        raise RuntimeError(
            "checkpoint_label_scalers missing from context. "
            "Run load_adme_checkpoint before transform_labels_from_checkpoint."
        )

    if transform_metadata is None:
        raise RuntimeError(
            "checkpoint_label_transform_metadata missing from context. "
            "Run load_adme_checkpoint before transform_labels_from_checkpoint."
        )

    missing_label_cols = [
        col for col in label_cols
        if col not in df.columns
    ]

    if missing_label_cols:
        raise KeyError(
            "The following label columns were requested but are missing from "
            f"the dataframe: {missing_label_cols}"
        )

    missing_scalers = [
        col for col in label_cols
        if col not in label_scalers
    ]

    if missing_scalers:
        raise KeyError(
            "The checkpoint does not contain label scalers for: "
            f"{missing_scalers}"
        )

    transformed_metadata_used = {}

    for col in label_cols:
        metadata = transform_metadata.get(col, {"transform": "identity"})
        scaler = label_scalers[col]

        raw_non_null = int(df[col].notna().sum())

        transformed = _apply_checkpoint_label_transform(df[col], metadata)

        transformed_non_null = int(transformed.notna().sum())

        df[col] = np.nan
        mask = transformed.notna()

        if mask.sum() > 0:
            scaled = scaler.transform(
                transformed.loc[mask].values.reshape(-1, 1)
            ).flatten()

            df.loc[mask, col] = scaled

        scaled_non_null = int(df[col].notna().sum())

        logger.debug(
            f"Checkpoint label transform '{col}' | "
            f"raw non-null: {raw_non_null}; "
            f"transformed non-null: {transformed_non_null}; "
            f"scaled non-null: {scaled_non_null}; "
            f"transform: {metadata.get('transform', 'identity')}"
        )

        transformed_metadata_used[col] = metadata

    context["dataframe"] = df
    context["task_names"] = label_cols
    context["label_scalers"] = label_scalers
    context["label_transform_metadata"] = transformed_metadata_used

    return context


def _is_finite_tensor(x):
    return torch.isfinite(x)


def _batch_to_device(batch, device):
    if hasattr(batch, "to"):
        return batch.to(device)

    return batch


def _model_forward_raw(model, batch):
    try:
        return model(batch)
    except TypeError:
        pass

    kwargs = {}

    for attr in [
        "x",
        "edge_index",
        "edge_attr",
        "batch",
        "global_features",
        "fp",
    ]:
        if hasattr(batch, attr):
            kwargs[attr] = getattr(batch, attr)

    try:
        return model(**kwargs)
    except TypeError:
        pass

    args = []

    for attr in [
        "x",
        "edge_index",
        "edge_attr",
        "batch",
        "global_features",
        "fp",
    ]:
        if hasattr(batch, attr):
            args.append(getattr(batch, attr))

    return model(*args)


def _extract_prediction_from_output(output):
    if torch.is_tensor(output):
        return output

    if isinstance(output, dict):
        for key in ["pred", "preds", "prediction", "predictions", "y_hat", "out", "output"]:
            value = output.get(key)
            if torch.is_tensor(value):
                return value

    if isinstance(output, (tuple, list)):
        for value in output:
            if torch.is_tensor(value):
                return value

    raise TypeError(
        "Could not extract prediction tensor from model output. "
        f"Output type was: {type(output)}"
    )


def _extract_embedding_from_output(output):
    if isinstance(output, dict):
        for key in ["embedding", "embeddings", "latent", "repr", "representation", "h", "z"]:
            value = output.get(key)
            if torch.is_tensor(value):
                return value

    if isinstance(output, (tuple, list)) and len(output) >= 2:
        for value in output[1:]:
            if torch.is_tensor(value):
                return value

    return None


def _predict(model, batch):
    output = _model_forward_raw(model, batch)
    return _extract_prediction_from_output(output)


def _predict_and_maybe_embed_auto(model, batch):
    try:
        output = model(batch, return_embeddings=True)
        pred = _extract_prediction_from_output(output)
        emb = _extract_embedding_from_output(output)
        return pred, emb
    except Exception:
        pass

    for method_name in [
        "forward_with_embeddings",
        "predict_with_embeddings",
        "get_predictions_and_embeddings",
    ]:
        if hasattr(model, method_name):
            method = getattr(model, method_name)
            try:
                output = method(batch)
                pred = _extract_prediction_from_output(output)
                emb = _extract_embedding_from_output(output)
                return pred, emb
            except Exception:
                pass

    pred_output = _model_forward_raw(model, batch)
    pred = _extract_prediction_from_output(pred_output)
    emb = _extract_embedding_from_output(pred_output)

    if emb is not None:
        return pred, emb

    for method_name in ["encode", "embed", "get_embeddings", "get_representation"]:
        if hasattr(model, method_name):
            method = getattr(model, method_name)
            try:
                emb = method(batch)
                if torch.is_tensor(emb):
                    return pred, emb
            except Exception:
                pass

    return pred, None


def _find_named_module(model, layer_name):
    if layer_name is None:
        return None

    for name, module in model.named_modules():
        if name == layer_name:
            return module

    return None


def _predict_and_embed_with_hook(model, batch, layer_name):
    module = _find_named_module(model, layer_name)

    if module is None:
        raise KeyError(
            f"Could not find embedding layer '{layer_name}' in model.named_modules()."
        )

    captured = {}

    def hook_fn(module, inputs, output):
        if torch.is_tensor(output):
            captured["embedding"] = output
        elif isinstance(output, (tuple, list)):
            for value in output:
                if torch.is_tensor(value):
                    captured["embedding"] = value
                    break
        elif isinstance(output, dict):
            for value in output.values():
                if torch.is_tensor(value):
                    captured["embedding"] = value
                    break

    handle = module.register_forward_hook(hook_fn)

    try:
        pred = _predict(model, batch)
    finally:
        handle.remove()

    emb = captured.get("embedding")

    return pred, emb


def _flatten_embedding_tensor(emb):
    if emb is None:
        return None

    if emb.dim() == 1:
        return emb.reshape(1, -1)

    if emb.dim() == 2:
        return emb

    if emb.dim() > 2:
        return emb.reshape(emb.shape[0], -1)

    return emb


def _get_batch_y(batch, num_tasks):
    if not hasattr(batch, "y"):
        raise AttributeError("Batch does not have attribute 'y'.")

    y = batch.y

    if y.dim() == 1:
        y = y.reshape(-1, num_tasks)

    elif y.dim() > 2:
        y = y.reshape(y.shape[0], -1)

    return y


def _get_batch_indices(batch):
    if hasattr(batch, "idx"):
        idx = batch.idx
        if torch.is_tensor(idx):
            return idx.detach().cpu().numpy().astype(int).tolist()
        return list(idx)

    return None


def _r2_score_np(y_true, y_pred):
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)

    if len(y_true) < 2:
        return np.nan

    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)

    if ss_tot <= 0:
        return np.nan

    return 1.0 - (ss_res / ss_tot)


def _pearson_np(y_true, y_pred):
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)

    if len(y_true) < 2:
        return np.nan

    if np.std(y_true) <= 0 or np.std(y_pred) <= 0:
        return np.nan

    return float(np.corrcoef(y_true, y_pred)[0, 1])


def _rmse_np(y_true, y_pred):
    return float(np.sqrt(np.mean((np.asarray(y_true) - np.asarray(y_pred)) ** 2)))


def _mae_np(y_true, y_pred):
    return float(np.mean(np.abs(np.asarray(y_true) - np.asarray(y_pred))))


def _collect_predictions_and_embeddings(
    model,
    loader,
    device,
    label_cols,
    embedding_cfg,
):
    model.eval()

    num_tasks = len(label_cols)

    all_pred = []
    all_y = []
    all_idx = []
    all_embeddings = []

    embedding_enabled = embedding_cfg.get("enabled", True)
    embedding_method = embedding_cfg.get("method", "auto")
    layer_name = embedding_cfg.get("layer_name")

    with torch.no_grad():
        for batch in loader:
            batch = _batch_to_device(batch, device)

            if embedding_enabled and embedding_method == "hook":
                pred, emb = _predict_and_embed_with_hook(
                    model=model,
                    batch=batch,
                    layer_name=layer_name,
                )
            elif embedding_enabled and layer_name is not None:
                pred, emb = _predict_and_embed_with_hook(
                    model=model,
                    batch=batch,
                    layer_name=layer_name,
                )
            elif embedding_enabled:
                pred, emb = _predict_and_maybe_embed_auto(model, batch)
            else:
                pred = _predict(model, batch)
                emb = None

            y = _get_batch_y(batch, num_tasks)

            if pred.dim() == 1:
                pred = pred.reshape(-1, num_tasks)

            all_pred.append(pred.detach().cpu())
            all_y.append(y.detach().cpu())

            idx = _get_batch_indices(batch)
            if idx is not None:
                all_idx.extend(idx)

            emb = _flatten_embedding_tensor(emb)

            if emb is not None:
                all_embeddings.append(emb.detach().cpu())

    pred = torch.cat(all_pred, dim=0).numpy()
    y = torch.cat(all_y, dim=0).numpy()

    if all_idx:
        idx = np.asarray(all_idx, dtype=int)
    else:
        idx = np.arange(pred.shape[0], dtype=int)

    if all_embeddings:
        embeddings = torch.cat(all_embeddings, dim=0).numpy()
    else:
        embeddings = None

    return {
        "pred": pred,
        "y": y,
        "idx": idx,
        "embeddings": embeddings,
    }


def _prediction_metrics_dataframe(split_name, collected, label_cols):
    pred = collected["pred"]
    y = collected["y"]

    rows = []

    for task_idx, task_name in enumerate(label_cols):
        y_task = y[:, task_idx]
        pred_task = pred[:, task_idx]

        mask = np.isfinite(y_task) & np.isfinite(pred_task)

        y_valid = y_task[mask]
        pred_valid = pred_task[mask]

        n = int(mask.sum())

        if n == 0:
            rows.append({
                "split": split_name,
                "task": task_name,
                "n": 0,
                "r2": np.nan,
                "pearson": np.nan,
                "mae": np.nan,
                "rmse": np.nan,
                "label_mean": np.nan,
                "label_std": np.nan,
                "pred_mean": np.nan,
                "pred_std": np.nan,
                "pred_std_over_label_std": np.nan,
                "collapse_flag_pred_std_lt_0_05_label_std": False,
            })
            continue

        label_std = float(np.std(y_valid))
        pred_std = float(np.std(pred_valid))

        ratio = (
            pred_std / label_std
            if label_std > 0 else np.nan
        )

        rows.append({
            "split": split_name,
            "task": task_name,
            "n": n,
            "r2": _r2_score_np(y_valid, pred_valid),
            "pearson": _pearson_np(y_valid, pred_valid),
            "mae": _mae_np(y_valid, pred_valid),
            "rmse": _rmse_np(y_valid, pred_valid),
            "label_mean": float(np.mean(y_valid)),
            "label_std": label_std,
            "pred_mean": float(np.mean(pred_valid)),
            "pred_std": pred_std,
            "pred_std_over_label_std": ratio,
            "collapse_flag_pred_std_lt_0_05_label_std": bool(
                np.isfinite(ratio) and ratio < 0.05
            ),
        })

    return pd.DataFrame(rows)


def _task_coverage_dataframe(df, label_cols, split_membership=None):
    rows = []

    for col in label_cols:
        row = {
            "task": col,
            "n_total_rows": len(df),
            "n_non_null": int(df[col].notna().sum()),
            "fraction_non_null": float(df[col].notna().mean()),
            "mean": _safe_float(df[col].mean(skipna=True)),
            "std": _safe_float(df[col].std(skipna=True)),
            "min": _safe_float(df[col].min(skipna=True)),
            "max": _safe_float(df[col].max(skipna=True)),
        }

        if split_membership is not None:
            for split_name, indices in split_membership.items():
                split_df = df.iloc[indices]
                row[f"{split_name}_n_rows"] = len(split_df)
                row[f"{split_name}_n_non_null"] = int(split_df[col].notna().sum())
                row[f"{split_name}_fraction_non_null"] = float(split_df[col].notna().mean())

        rows.append(row)

    return pd.DataFrame(rows)


def _task_overlap_matrix(df, label_cols):
    overlap = pd.DataFrame(
        0,
        index=label_cols,
        columns=label_cols,
        dtype=int,
    )

    for col_i in label_cols:
        mask_i = df[col_i].notna()

        for col_j in label_cols:
            mask_j = df[col_j].notna()
            overlap.loc[col_i, col_j] = int((mask_i & mask_j).sum())

    return overlap


def _task_overlap_fraction_matrix(df, label_cols):
    overlap_fraction = pd.DataFrame(
        np.nan,
        index=label_cols,
        columns=label_cols,
        dtype=float,
    )

    for col_i in label_cols:
        mask_i = df[col_i].notna()

        for col_j in label_cols:
            mask_j = df[col_j].notna()
            denom = int(mask_i.sum())

            if denom > 0:
                overlap_fraction.loc[col_i, col_j] = float((mask_i & mask_j).sum() / denom)

    return overlap_fraction


def _task_group_summary(metrics_df, task_groups):
    if not task_groups:
        return pd.DataFrame()

    rows = []

    for split_name in sorted(metrics_df["split"].unique()):
        split_df = metrics_df[metrics_df["split"] == split_name]

        for group_name, tasks in task_groups.items():
            group_df = split_df[split_df["task"].isin(tasks)]

            rows.append({
                "split": split_name,
                "group": group_name,
                "n_tasks_configured": len(tasks),
                "n_tasks_present": int(len(group_df)),
                "median_r2": _safe_float(group_df["r2"].median(skipna=True)),
                "mean_r2": _safe_float(group_df["r2"].mean(skipna=True)),
                "median_pearson": _safe_float(group_df["pearson"].median(skipna=True)),
                "median_pred_std_over_label_std": _safe_float(
                    group_df["pred_std_over_label_std"].median(skipna=True)
                ),
                "n_collapse_flags": int(
                    group_df["collapse_flag_pred_std_lt_0_05_label_std"].sum()
                ),
                "total_labelled_points": int(group_df["n"].sum()),
            })

    return pd.DataFrame(rows)


def _effective_rank_from_embeddings(embeddings):
    if embeddings is None:
        return {
            "available": False,
            "reason": "No embeddings extracted.",
        }

    if embeddings.shape[0] < 2:
        return {
            "available": False,
            "reason": "Fewer than two embedding rows.",
        }

    x = embeddings.astype(float)
    x = x - np.nanmean(x, axis=0, keepdims=True)
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)

    try:
        singular_values = np.linalg.svd(x, compute_uv=False)
    except np.linalg.LinAlgError as exc:
        return {
            "available": False,
            "reason": f"SVD failed: {exc}",
        }

    singular_values = singular_values[singular_values > 0]

    if len(singular_values) == 0:
        return {
            "available": True,
            "n_rows": int(embeddings.shape[0]),
            "embedding_dim": int(embeddings.shape[1]),
            "effective_rank": 0.0,
            "participation_ratio": 0.0,
            "top_singular_value_fraction": np.nan,
        }

    p = singular_values / singular_values.sum()
    entropy = -np.sum(p * np.log(p + 1e-12))
    effective_rank = float(np.exp(entropy))

    participation_ratio = float(
        (np.sum(singular_values) ** 2) /
        np.sum(singular_values ** 2)
    )

    return {
        "available": True,
        "n_rows": int(embeddings.shape[0]),
        "embedding_dim": int(embeddings.shape[1]),
        "effective_rank": effective_rank,
        "participation_ratio": participation_ratio,
        "top_singular_value_fraction": float(singular_values[0] / singular_values.sum()),
        "n_nonzero_singular_values": int(len(singular_values)),
    }


def _embedding_pca_dataframe(embeddings, indices, n_components=10, max_rows=50000):
    if embeddings is None:
        return pd.DataFrame()

    if embeddings.shape[0] < 2:
        return pd.DataFrame()

    x = embeddings.astype(float)
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)

    if x.shape[0] > max_rows:
        rng = np.random.default_rng(42)
        keep = np.sort(rng.choice(np.arange(x.shape[0]), size=max_rows, replace=False))
        x = x[keep]
        indices = indices[keep]

    n_components = min(n_components, x.shape[0], x.shape[1])

    if n_components < 1:
        return pd.DataFrame()

    pca = PCA(n_components=n_components)
    pcs = pca.fit_transform(x)

    rows = []

    for i in range(pcs.shape[0]):
        row = {
            "row_index": int(indices[i]),
        }

        for j in range(pcs.shape[1]):
            row[f"PC{j + 1}"] = float(pcs[i, j])

        rows.append(row)

    return pd.DataFrame(rows)


def _parameter_norms_dataframe(model):
    rows = []

    for name, param in model.named_parameters():
        if param is None:
            continue

        data = param.detach().cpu()

        rows.append({
            "parameter": name,
            "shape": str(tuple(data.shape)),
            "numel": int(data.numel()),
            "norm_l2": float(torch.norm(data).item()),
            "mean": float(data.float().mean().item()),
            "std": float(data.float().std().item()) if data.numel() > 1 else 0.0,
            "abs_max": float(data.abs().max().item()) if data.numel() > 0 else 0.0,
        })

    return pd.DataFrame(rows)


def _task_head_norms_dataframe(model, label_cols):
    rows = []
    num_tasks = len(label_cols)

    for module_name, module in model.named_modules():
        if not isinstance(module, torch.nn.Linear):
            continue

        weight = module.weight.detach().cpu()
        bias = module.bias.detach().cpu() if module.bias is not None else None

        if weight.shape[0] == num_tasks:
            for task_idx, task_name in enumerate(label_cols):
                rows.append({
                    "module": module_name,
                    "task": task_name,
                    "orientation": "out_features_match_num_tasks",
                    "weight_norm_l2": float(torch.norm(weight[task_idx]).item()),
                    "weight_mean": float(weight[task_idx].float().mean().item()),
                    "weight_std": float(weight[task_idx].float().std().item())
                    if weight[task_idx].numel() > 1 else 0.0,
                    "bias": float(bias[task_idx].item()) if bias is not None else np.nan,
                })

        elif weight.shape[1] == num_tasks:
            for task_idx, task_name in enumerate(label_cols):
                rows.append({
                    "module": module_name,
                    "task": task_name,
                    "orientation": "in_features_match_num_tasks",
                    "weight_norm_l2": float(torch.norm(weight[:, task_idx]).item()),
                    "weight_mean": float(weight[:, task_idx].float().mean().item()),
                    "weight_std": float(weight[:, task_idx].float().std().item())
                    if weight[:, task_idx].numel() > 1 else 0.0,
                    "bias": np.nan,
                })

    return pd.DataFrame(rows)

def _task_groups_names_to_indices(task_groups, label_cols):
    if task_groups is None:
        return None

    label_to_idx = {
        task_name: task_idx
        for task_idx, task_name in enumerate(label_cols)
    }

    converted = {}

    for group_name, tasks in task_groups.items():
        converted_tasks = []

        for task in tasks:
            if isinstance(task, int):
                task_idx = task

            elif isinstance(task, str):
                if task not in label_to_idx:
                    raise KeyError(
                        f"Task '{task}' in task group '{group_name}' is not present "
                        f"in label_cols. Available tasks: {label_cols}"
                    )

                task_idx = label_to_idx[task]

            else:
                raise TypeError(
                    f"Task entries in task_groups must be strings or integers. "
                    f"Got {type(task)} in group '{group_name}'."
                )

            converted_tasks.append(task_idx)

        converted[group_name] = converted_tasks

    flattened = [
        task_idx
        for tasks in converted.values()
        for task_idx in tasks
    ]

    expected = list(range(len(label_cols)))
    got = sorted(flattened)

    if got != expected:
        missing = sorted(set(expected) - set(flattened))
        duplicated = sorted({
            task_idx
            for task_idx in flattened
            if flattened.count(task_idx) > 1
        })
        extra = sorted(set(flattened) - set(expected))

        raise ValueError(
            "Converted task_groups do not cover label_cols exactly once.\n"
            f"Missing task indices: {missing}\n"
            f"Duplicated task indices: {duplicated}\n"
            f"Extra task indices: {extra}\n"
            f"Expected indices: {expected}\n"
            f"Got indices: {got}"
        )

    return converted

def _infer_group_trunk_dims_from_state_dict(state_dict):
    """
    Infer per-task-group trunk dimensions from checkpoint tensor shapes.

    Expected checkpoint keys look like:
        group_trunks.physchem.0.weight
        group_trunks.physchem.3.weight

    For a Linear layer, weight shape is:
        [out_features, in_features]

    We use group_trunks.<group>.0.weight out_features as the group's trunk width.
    """

    group_trunk_dims = {}

    prefix = "group_trunks."

    for key, tensor in state_dict.items():
        if not key.startswith(prefix):
            continue

        if not key.endswith(".0.weight"):
            continue

        if not torch.is_tensor(tensor):
            continue

        parts = key.split(".")

        if len(parts) < 4:
            continue

        group_name = parts[1]
        out_features = int(tensor.shape[0])

        group_trunk_dims[group_name] = out_features

    return group_trunk_dims


def _infer_head_input_dims_from_state_dict(state_dict):
    """
    Infer task-head input dimensions from checkpoint tensor shapes.

    Expected checkpoint keys look like:
        heads.14.0.weight

    For a Linear layer, weight shape is:
        [out_features, in_features]

    The in_features dimension is the incoming group trunk width for that task.
    """

    head_input_dims = {}

    prefix = "heads."

    for key, tensor in state_dict.items():
        if not key.startswith(prefix):
            continue

        if not key.endswith(".0.weight"):
            continue

        if not torch.is_tensor(tensor):
            continue

        parts = key.split(".")

        if len(parts) < 4:
            continue

        try:
            task_idx = int(parts[1])
        except ValueError:
            continue

        in_features = int(tensor.shape[1])
        head_input_dims[task_idx] = in_features

    return head_input_dims


def _add_if_constructor_accepts(model_kwargs, accepted_args, candidate_names, value):
    """
    Add value to model_kwargs using the first candidate name accepted by the
    model constructor.
    """

    if value is None:
        return False

    for name in candidate_names:
        if name in accepted_args and name not in model_kwargs:
            model_kwargs[name] = value
            return True

    return False

def _group_architecture_candidates_from_dims(group_dims):
    """
    Build plausible group_architecture formats for older/newer model variants.

    The checkpoint tells us each group trunk output dim, e.g.
        physchem -> 128
        metabolism -> 512

    Different versions of GINRegressor may expect this as:
        {"physchem": 128}
    or:
        {"physchem": {"hidden_dim": 128, "output_dim": 128}}
    or:
        {"physchem": (128, 128)}
    """

    if not group_dims:
        return []

    candidates = []

    # Most likely if constructor uses group_architecture directly.
    candidates.append(group_dims)

    # Common explicit dict format.
    candidates.append({
        group: {
            "hidden_dim": dim,
            "output_dim": dim,
        }
        for group, dim in group_dims.items()
    })

    # Alternative explicit dict format.
    candidates.append({
        group: {
            "trunk_hidden_dim": dim,
            "trunk_output_dim": dim,
        }
        for group, dim in group_dims.items()
    })

    # Tuple format, if older code used unpacking.
    candidates.append({
        group: (dim, dim)
        for group, dim in group_dims.items()
    })

    return candidates

def _instantiate_model_from_checkpoint(config, context):
    ModelClass = context["model_class"]

    checkpoint = context["adme_checkpoint"]
    state_dict = _strip_module_prefix_from_state_dict(
        _get_state_dict_from_checkpoint(checkpoint)
    )

    label_cols = context.get("task_names")

    if label_cols is None:
        raise RuntimeError(
            "context['task_names'] is missing. Run transform_labels_from_checkpoint "
            "before diagnose_adme_checkpoint."
        )

    # Start with YAML model kwargs.
    model_kwargs = dict(config.get("model_kwargs", {}))

    # If the checkpoint contains saved model kwargs/hparams, use them as defaults.
    # YAML still wins if the same key is explicitly supplied there.
    checkpoint_model_kwargs = _checkpoint_get(
        checkpoint,
        keys=[
            "model_kwargs",
            "model_config",
            "model_hparams",
            "hyperparameters",
            "hparams",
            "best_params",
        ],
        default=None,
    )

    if isinstance(checkpoint_model_kwargs, dict):
        for key, value in checkpoint_model_kwargs.items():
            model_kwargs.setdefault(key, value)

    inferred_kwargs = {
        "input_dim": context.get("input_dim"),
        "edge_dim": context.get("edge_dim"),
        "global_feat_dim": context.get("global_feat_dim"),
        "fp_dim": context.get("fp_dim"),
        "num_tasks": context.get("num_tasks"),
    }

    signature = inspect.signature(ModelClass.__init__)
    accepted_args = set(signature.parameters.keys())
    accepted_args.discard("self")

    logger.info(f"Model constructor accepted args: {sorted(accepted_args)}")

    for key, value in inferred_kwargs.items():
        if key in accepted_args and key not in model_kwargs and value is not None:
            model_kwargs[key] = value

    if "task_groups" in accepted_args:
        raw_task_groups = model_kwargs.get("task_groups")

        if raw_task_groups is None:
            raw_task_groups = config.get("task_groups")

        if raw_task_groups is None:
            raise ValueError(
                "Model constructor requires task_groups, but none were provided. "
                "Add diagnose_adme_checkpoint.task_groups to the YAML."
            )

        model_kwargs["task_groups"] = _task_groups_names_to_indices(
            raw_task_groups,
            label_cols,
        )

    if (
        "task_names" in accepted_args
        and "task_names" not in model_kwargs
    ):
        model_kwargs["task_names"] = label_cols

    if (
        "label_cols" in accepted_args
        and "label_cols" not in model_kwargs
    ):
        model_kwargs["label_cols"] = label_cols

    # Infer group-specific trunk widths from checkpoint.
    inferred_group_trunk_dims = _infer_group_trunk_dims_from_state_dict(state_dict)
    inferred_head_input_dims = _infer_head_input_dims_from_state_dict(state_dict)

    logger.info(f"Inferred group trunk dims from checkpoint: {inferred_group_trunk_dims}")
    logger.info(f"Inferred head input dims from checkpoint: {inferred_head_input_dims}")

    # Try common constructor names for per-group hidden/output dimensions.
    added_group_dims = _add_if_constructor_accepts(
        model_kwargs=model_kwargs,
        accepted_args=accepted_args,
        candidate_names=[
            "group_trunk_dims",
            "group_hidden_dims",
            "group_output_dims",
            "task_group_dims",
            "task_group_hidden_dims",
            "task_group_output_dims",
            "group_dims",
        ],
        value=inferred_group_trunk_dims,
    )

    if not added_group_dims and inferred_group_trunk_dims:
        logger.warning(
            "Checkpoint contains group-specific trunk dimensions, but the model "
            "constructor does not appear to accept any known group-dimension kwarg. "
            "Known candidates tried: group_trunk_dims, group_hidden_dims, "
            "group_output_dims, task_group_dims, task_group_hidden_dims, "
            "task_group_output_dims, group_dims."
        )

    strict_load = config.get("strict_load", False)

    architecture_attempts = [None]

    if (
        "group_architecture" in accepted_args
        and "group_architecture" not in model_kwargs
        and inferred_group_trunk_dims
    ):
        architecture_attempts = _group_architecture_candidates_from_dims(
            inferred_group_trunk_dims
        )

    last_error = None

    for attempt_idx, group_architecture in enumerate(architecture_attempts):

        trial_kwargs = dict(model_kwargs)

        if group_architecture is not None:
            trial_kwargs["group_architecture"] = group_architecture

        logger.info("=" * 80)
        logger.info(
            f"Trying model instantiation attempt {attempt_idx + 1}/"
            f"{len(architecture_attempts)}"
        )
        logger.info(f"Instantiating model class: {ModelClass}")
        logger.info(f"Model kwargs: {trial_kwargs}")

        try:
            model = ModelClass(**trial_kwargs)

            load_result = model.load_state_dict(
                state_dict,
                strict=strict_load,
            )

            if hasattr(load_result, "missing_keys"):
                logger.info(f"Checkpoint missing keys: {load_result.missing_keys}")
                logger.info(f"Checkpoint unexpected keys: {load_result.unexpected_keys}")

            logger.info("Successfully instantiated model and loaded checkpoint.")
            return model

        except Exception as exc:
            last_error = exc

            logger.warning(
                f"Model instantiation/load attempt {attempt_idx + 1} failed: {exc}"
            )

    logger.error("Checkpoint loading failed for all architecture attempts.")
    logger.error(
        "This means the diagnostic model kwargs still do not exactly match the "
        "training-time model architecture."
    )
    logger.error(f"Inferred group trunk dims: {inferred_group_trunk_dims}")
    logger.error(f"Inferred head input dims: {inferred_head_input_dims}")
    logger.error(f"Base model kwargs used: {model_kwargs}")

    raise last_error

def _shared_parameter_vector_from_grads(model, exclude_patterns):
    chunks = []

    for name, param in model.named_parameters():
        if param.grad is None:
            continue

        lowered = name.lower()

        if any(pattern.lower() in lowered for pattern in exclude_patterns):
            continue

        chunks.append(param.grad.detach().flatten().cpu())

    if not chunks:
        return None

    return torch.cat(chunks, dim=0)


def _gradient_cosine_matrix(
    model,
    loader,
    device,
    label_cols,
    max_batches=4,
    exclude_patterns=None,
):
    if exclude_patterns is None:
        exclude_patterns = ["head", "heads", "task", "output", "readout"]

    num_tasks = len(label_cols)

    grad_sums = {
        task_idx: None
        for task_idx in range(num_tasks)
    }

    task_counts = {
        task_idx: 0
        for task_idx in range(num_tasks)
    }

    model.train()

    batches_used = 0

    for batch in loader:
        if batches_used >= max_batches:
            break

        batch = _batch_to_device(batch, device)
        y = _get_batch_y(batch, num_tasks)

        for task_idx in range(num_tasks):
            mask = torch.isfinite(y[:, task_idx])

            if int(mask.sum().item()) == 0:
                continue

            model.zero_grad(set_to_none=True)

            pred = _predict(model, batch)

            if pred.dim() == 1:
                pred = pred.reshape(-1, num_tasks)

            loss = torch.nn.functional.mse_loss(
                pred[mask, task_idx],
                y[mask, task_idx],
            )

            loss.backward()

            grad_vec = _shared_parameter_vector_from_grads(
                model,
                exclude_patterns=exclude_patterns,
            )

            if grad_vec is None:
                continue

            if grad_sums[task_idx] is None:
                grad_sums[task_idx] = grad_vec
            else:
                min_len = min(len(grad_sums[task_idx]), len(grad_vec))
                grad_sums[task_idx] = grad_sums[task_idx][:min_len] + grad_vec[:min_len]

            task_counts[task_idx] += 1

        batches_used += 1

    rows = []

    for i, task_i in enumerate(label_cols):
        row = {
            "task": task_i,
        }

        grad_i = grad_sums[i]

        for j, task_j in enumerate(label_cols):
            grad_j = grad_sums[j]

            if grad_i is None or grad_j is None:
                cosine = np.nan
            else:
                min_len = min(len(grad_i), len(grad_j))
                gi = grad_i[:min_len]
                gj = grad_j[:min_len]

                denom = torch.norm(gi) * torch.norm(gj)

                if denom.item() <= 0:
                    cosine = np.nan
                else:
                    cosine = float(torch.dot(gi, gj).item() / denom.item())

            row[task_j] = cosine

        rows.append(row)

    matrix_df = pd.DataFrame(rows)

    count_df = pd.DataFrame([
        {
            "task": label_cols[task_idx],
            "n_batches_with_labels": task_counts[task_idx],
        }
        for task_idx in range(num_tasks)
    ])

    model.eval()

    return matrix_df, count_df


def _split_indices_from_graphs(graph_list):
    indices = []

    for graph in graph_list:
        if hasattr(graph, "idx"):
            idx = graph.idx

            if torch.is_tensor(idx):
                idx = int(idx.item())
            else:
                idx = int(idx)

            indices.append(idx)

    return indices


@register_task(
    "diagnose_adme_checkpoint",
    category="ADME",
    description="Run checkpoint diagnostics: task coverage, predictions, embeddings, head norms and gradient conflict."
)
def diagnose_adme_checkpoint(config, context):
    out_dir = _ensure_dir(config.get("out_dir", "outputs/mtl_adme/diagnostics"))
    device = torch.device(config.get("device", context.get("device", "cpu")))

    label_cols = context.get("task_names")
    if label_cols is None:
        raise RuntimeError(
            "context['task_names'] is missing. Run transform_labels_from_checkpoint first."
        )

    train_loader = context.get("train_loader")
    val_loader = context.get("val_loader")

    if train_loader is None or val_loader is None:
        raise RuntimeError(
            "train_loader and val_loader are required. Run split_data before diagnostics."
        )

    model = _instantiate_model_from_checkpoint(config, context)
    model = model.to(device)
    model.eval()

    df = context["dataframe"]

    train_indices = _split_indices_from_graphs(context.get("train_list", []))
    val_indices = _split_indices_from_graphs(context.get("val_list", []))

    split_membership = {}

    if train_indices:
        split_membership["train"] = train_indices

    if val_indices:
        split_membership["val"] = val_indices

    logger.info("=" * 80)
    logger.info("Running ADME checkpoint diagnostics")
    logger.info(f"Output directory: {out_dir}")
    logger.info(f"Number of tasks: {len(label_cols)}")

    coverage_df = _task_coverage_dataframe(
        df=df,
        label_cols=label_cols,
        split_membership=split_membership,
    )

    coverage_path = os.path.join(out_dir, "task_label_coverage.csv")
    coverage_df.to_csv(coverage_path, index=False)
    logger.info(f"Saved task label coverage: {coverage_path}")

    overlap_counts = _task_overlap_matrix(df, label_cols)
    overlap_counts_path = os.path.join(out_dir, "task_overlap_counts.csv")
    overlap_counts.to_csv(overlap_counts_path)
    logger.info(f"Saved task overlap count matrix: {overlap_counts_path}")

    overlap_fraction = _task_overlap_fraction_matrix(df, label_cols)
    overlap_fraction_path = os.path.join(out_dir, "task_overlap_fraction_by_row_task.csv")
    overlap_fraction.to_csv(overlap_fraction_path)
    logger.info(f"Saved task overlap fraction matrix: {overlap_fraction_path}")

    parameter_norms_df = _parameter_norms_dataframe(model)
    parameter_norms_path = os.path.join(out_dir, "model_parameter_norms.csv")
    parameter_norms_df.to_csv(parameter_norms_path, index=False)
    logger.info(f"Saved model parameter norms: {parameter_norms_path}")

    task_head_norms_df = _task_head_norms_dataframe(model, label_cols)
    task_head_norms_path = os.path.join(out_dir, "task_head_norms.csv")
    task_head_norms_df.to_csv(task_head_norms_path, index=False)
    logger.info(f"Saved task head norms: {task_head_norms_path}")

    prediction_splits = config.get("prediction_splits", ["train", "val"])
    embedding_cfg = config.get("embedding", {"enabled": True, "method": "auto"})

    loader_map = {
        "train": train_loader,
        "val": val_loader,
    }

    all_metrics = []
    embedding_rank_summary = {}

    for split_name in prediction_splits:
        if split_name not in loader_map:
            raise KeyError(
                f"Unknown prediction split '{split_name}'. Supported: train, val."
            )

        logger.info(f"Collecting predictions for split: {split_name}")

        collected = _collect_predictions_and_embeddings(
            model=model,
            loader=loader_map[split_name],
            device=device,
            label_cols=label_cols,
            embedding_cfg=embedding_cfg,
        )

        pred_df = pd.DataFrame(
            collected["pred"],
            columns=[f"pred_{task}" for task in label_cols],
        )
        pred_df.insert(0, "row_index", collected["idx"])

        y_df = pd.DataFrame(
            collected["y"],
            columns=[f"label_{task}" for task in label_cols],
        )
        y_df.insert(0, "row_index", collected["idx"])

        pred_path = os.path.join(out_dir, f"{split_name}_predictions_scaled.csv")
        labels_path = os.path.join(out_dir, f"{split_name}_labels_scaled.csv")

        pred_df.to_csv(pred_path, index=False)
        y_df.to_csv(labels_path, index=False)

        logger.info(f"Saved {split_name} predictions: {pred_path}")
        logger.info(f"Saved {split_name} labels: {labels_path}")

        metrics_df = _prediction_metrics_dataframe(
            split_name=split_name,
            collected=collected,
            label_cols=label_cols,
        )

        all_metrics.append(metrics_df)

        embeddings = collected["embeddings"]

        embedding_rank_summary[split_name] = _effective_rank_from_embeddings(embeddings)

        if embeddings is not None:
            pca_components = embedding_cfg.get("pca_components", 10)
            max_rows_for_pca = embedding_cfg.get("max_rows_for_pca", 50000)

            pca_df = _embedding_pca_dataframe(
                embeddings=embeddings,
                indices=collected["idx"],
                n_components=pca_components,
                max_rows=max_rows_for_pca,
            )

            pca_path = os.path.join(out_dir, f"{split_name}_embedding_pca.csv")
            pca_df.to_csv(pca_path, index=False)
            logger.info(f"Saved {split_name} embedding PCA: {pca_path}")

        else:
            logger.warning(
                f"No embeddings were extracted for split '{split_name}'. "
                "If your model does not expose embeddings automatically, set "
                "diagnose_adme_checkpoint.embedding.layer_name to a named module."
            )

    metrics_all_df = pd.concat(all_metrics, axis=0, ignore_index=True)

    metrics_path = os.path.join(out_dir, "task_prediction_metrics.csv")
    metrics_all_df.to_csv(metrics_path, index=False)
    logger.info(f"Saved task prediction metrics: {metrics_path}")

    task_groups = config.get("task_groups", {})
    group_summary_df = _task_group_summary(metrics_all_df, task_groups)

    group_summary_path = os.path.join(out_dir, "task_group_prediction_summary.csv")
    group_summary_df.to_csv(group_summary_path, index=False)
    logger.info(f"Saved task group prediction summary: {group_summary_path}")

    embedding_rank_path = os.path.join(out_dir, "embedding_effective_rank.json")
    with open(embedding_rank_path, "w") as f:
        json.dump(embedding_rank_summary, f, indent=2)

    logger.info(f"Saved embedding effective-rank summary: {embedding_rank_path}")

    gradient_cfg = config.get("gradients", {})
    gradient_enabled = gradient_cfg.get("enabled", True)

    if gradient_enabled:
        gradient_split = gradient_cfg.get("split", "train")
        max_batches = gradient_cfg.get("max_batches", 4)
        exclude_patterns = gradient_cfg.get(
            "exclude_param_name_patterns",
            ["head", "heads", "task", "output", "readout"],
        )

        if gradient_split not in loader_map:
            raise KeyError(
                f"Unknown gradient split '{gradient_split}'. Supported: train, val."
            )

        logger.info(
            f"Computing task gradient cosine matrix on split '{gradient_split}' "
            f"using up to {max_batches} batches"
        )

        grad_cosine_df, grad_count_df = _gradient_cosine_matrix(
            model=model,
            loader=loader_map[gradient_split],
            device=device,
            label_cols=label_cols,
            max_batches=max_batches,
            exclude_patterns=exclude_patterns,
        )

        grad_cosine_path = os.path.join(out_dir, "task_gradient_cosine_matrix.csv")
        grad_count_path = os.path.join(out_dir, "task_gradient_batch_counts.csv")

        grad_cosine_df.to_csv(grad_cosine_path, index=False)
        grad_count_df.to_csv(grad_count_path, index=False)

        logger.info(f"Saved task gradient cosine matrix: {grad_cosine_path}")
        logger.info(f"Saved task gradient batch counts: {grad_count_path}")

    summary = {
        "out_dir": out_dir,
        "task_label_coverage": coverage_path,
        "task_overlap_counts": overlap_counts_path,
        "task_overlap_fraction": overlap_fraction_path,
        "task_prediction_metrics": metrics_path,
        "task_group_prediction_summary": group_summary_path,
        "model_parameter_norms": parameter_norms_path,
        "task_head_norms": task_head_norms_path,
        "embedding_effective_rank": embedding_rank_path,
    }

    summary_path = os.path.join(out_dir, "diagnostic_outputs.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    logger.info(f"Saved diagnostic output manifest: {summary_path}")
    logger.info("Finished ADME checkpoint diagnostics")
    logger.info("=" * 80)

    context["diagnostic_summary"] = summary

    return context

@register_task(
    "plot_adme_checkpoint_diagnostics",
    category="ADME",
    description="Generate diagnostic plots from ADME checkpoint diagnostic CSV outputs."
)
def plot_adme_checkpoint_diagnostics(config, context):
    diag_dir = config.get(
        "diag_dir",
        "outputs/mtl_adme/diagnostics/checkpoint_diagnostics",
    )

    plot_dir = config.get(
        "plot_dir",
        os.path.join(diag_dir, "plots"),
    )

    os.makedirs(plot_dir, exist_ok=True)

    split_for_task_bars = config.get("split_for_task_bars", "val")
    scatter_split = config.get("scatter_split", "val")

    worst_n_scatter_tasks = int(config.get("worst_n_scatter_tasks", 12))
    best_n_scatter_tasks = int(config.get("best_n_scatter_tasks", 6))

    heatmap_figsize = tuple(config.get("heatmap_figsize", [14, 12]))
    bar_figsize = tuple(config.get("bar_figsize", [14, 8]))
    scatter_figsize = tuple(config.get("scatter_figsize", [16, 12]))
    dpi = int(config.get("dpi", 200))

    logger.info("=" * 80)
    logger.info("Generating ADME checkpoint diagnostic plots")
    logger.info(f"Diagnostic CSV directory: {diag_dir}")
    logger.info(f"Plot directory: {plot_dir}")

    def path(name):
        return os.path.join(diag_dir, name)

    def savefig(filename):
        out_path = os.path.join(plot_dir, filename)
        plt.tight_layout()
        plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved plot: {out_path}")
        return out_path

    def wrap_label(label, width=28):
        return "\n".join(textwrap.wrap(str(label), width=width))

    def read_csv_if_exists(filename, **kwargs):
        file_path = path(filename)

        if not os.path.exists(file_path):
            logger.warning(f"Missing diagnostic file, skipping: {file_path}")
            return None

        try:
            return pd.read_csv(file_path, **kwargs)
        except pd.errors.EmptyDataError:
            logger.warning(f"Diagnostic file is empty, skipping: {file_path}")
            return None

    outputs = {}

    metrics_df = read_csv_if_exists("task_prediction_metrics.csv")
    coverage_df = read_csv_if_exists("task_label_coverage.csv")
    group_summary_df = read_csv_if_exists("task_group_prediction_summary.csv")
    overlap_counts_df = read_csv_if_exists("task_overlap_counts.csv", index_col=0)
    overlap_fraction_df = read_csv_if_exists(
        "task_overlap_fraction_by_row_task.csv",
        index_col=0,
    )
    grad_cosine_df = read_csv_if_exists("task_gradient_cosine_matrix.csv")
    grad_counts_df = read_csv_if_exists("task_gradient_batch_counts.csv")

    # -------------------------------------------------------------------------
    # 1. Validation R2 by task
    # -------------------------------------------------------------------------
    if metrics_df is not None:
        split_df = metrics_df[metrics_df["split"] == split_for_task_bars].copy()

        if len(split_df) > 0:
            split_df = split_df.sort_values("r2", ascending=True)

            plt.figure(figsize=bar_figsize)

            colors = [
                "tab:red" if pd.notna(v) and v < 0 else
                "tab:orange" if pd.notna(v) and v < 0.2 else
                "tab:blue"
                for v in split_df["r2"]
            ]

            plt.barh(
                [wrap_label(t) for t in split_df["task"]],
                split_df["r2"],
                color=colors,
            )

            plt.axvline(0.0, color="black", linewidth=1)
            plt.axvline(0.2, color="grey", linewidth=1, linestyle="--")
            plt.xlabel("R²")
            plt.ylabel("Task")
            plt.title(f"{split_for_task_bars} R² by task")
            outputs["r2_by_task"] = savefig(
                f"{split_for_task_bars}_r2_by_task.png"
            )

    # -------------------------------------------------------------------------
    # 2. Prediction std / label std by task
    # -------------------------------------------------------------------------
    if metrics_df is not None:
        split_df = metrics_df[metrics_df["split"] == split_for_task_bars].copy()

        if len(split_df) > 0:
            split_df = split_df.sort_values(
                "pred_std_over_label_std",
                ascending=True,
            )

            plt.figure(figsize=bar_figsize)

            ratio = split_df["pred_std_over_label_std"].replace(
                [np.inf, -np.inf],
                np.nan,
            )

            colors = [
                "tab:red" if pd.notna(v) and v < 0.05 else
                "tab:orange" if pd.notna(v) and v < 0.25 else
                "tab:blue"
                for v in ratio
            ]

            plt.barh(
                [wrap_label(t) for t in split_df["task"]],
                ratio,
                color=colors,
            )

            plt.axvline(0.05, color="red", linewidth=1, linestyle="--")
            plt.axvline(0.25, color="orange", linewidth=1, linestyle="--")
            plt.axvline(1.0, color="black", linewidth=1)
            plt.xlabel("Prediction standard deviation / label standard deviation")
            plt.ylabel("Task")
            plt.title(
                f"{split_for_task_bars} prediction variance ratio by task\n"
                "Low values indicate mean-like or collapsed predictions"
            )
            outputs["prediction_variance_ratio"] = savefig(
                f"{split_for_task_bars}_prediction_std_over_label_std.png"
            )

    # -------------------------------------------------------------------------
    # 3. Train vs validation R2 comparison
    # -------------------------------------------------------------------------
    if metrics_df is not None:
        train_df = metrics_df[metrics_df["split"] == "train"][["task", "r2"]]
        val_df = metrics_df[metrics_df["split"] == "val"][["task", "r2"]]

        if len(train_df) > 0 and len(val_df) > 0:
            merged = train_df.merge(
                val_df,
                on="task",
                suffixes=("_train", "_val"),
            )

            plt.figure(figsize=(8, 8))
            plt.scatter(
                merged["r2_train"],
                merged["r2_val"],
                alpha=0.8,
            )

            min_val = np.nanmin([
                merged["r2_train"].min(),
                merged["r2_val"].min(),
                0.0,
            ])
            max_val = np.nanmax([
                merged["r2_train"].max(),
                merged["r2_val"].max(),
                1.0,
            ])

            plt.plot([min_val, max_val], [min_val, max_val], color="black")
            plt.axhline(0.0, color="grey", linestyle="--", linewidth=1)
            plt.axvline(0.0, color="grey", linestyle="--", linewidth=1)
            plt.xlabel("Train R²")
            plt.ylabel("Validation R²")
            plt.title("Train vs validation R² by task")

            for _, row in merged.iterrows():
                if (
                    pd.notna(row["r2_val"]) and row["r2_val"] < 0.1
                ) or (
                    pd.notna(row["r2_train"]) and
                    pd.notna(row["r2_val"]) and
                    row["r2_train"] - row["r2_val"] > 0.5
                ):
                    plt.annotate(
                        row["task"],
                        (row["r2_train"], row["r2_val"]),
                        fontsize=6,
                        alpha=0.8,
                    )

            outputs["train_vs_val_r2"] = savefig("train_vs_val_r2.png")

    # -------------------------------------------------------------------------
    # 4. Task label coverage
    # -------------------------------------------------------------------------
    if coverage_df is not None:
        coverage_plot_df = coverage_df.sort_values("n_non_null", ascending=True)

        plt.figure(figsize=bar_figsize)
        plt.barh(
            [wrap_label(t) for t in coverage_plot_df["task"]],
            coverage_plot_df["n_non_null"],
            color="tab:blue",
        )
        plt.xlabel("Number of labelled compounds")
        plt.ylabel("Task")
        plt.title("Task label coverage after filtering/transforms")
        outputs["task_label_coverage"] = savefig("task_label_coverage.png")

    # -------------------------------------------------------------------------
    # 5. Group-level summary
    # -------------------------------------------------------------------------
    if group_summary_df is not None and len(group_summary_df) > 0:
        group_split_df = group_summary_df[
            group_summary_df["split"] == split_for_task_bars
        ].copy()

        if len(group_split_df) > 0:
            group_split_df = group_split_df.sort_values("median_r2")

            plt.figure(figsize=(12, 6))
            plt.bar(
                group_split_df["group"],
                group_split_df["median_r2"],
                color="tab:purple",
            )
            plt.axhline(0.0, color="black", linewidth=1)
            plt.axhline(0.2, color="grey", linestyle="--", linewidth=1)
            plt.xticks(rotation=45, ha="right")
            plt.ylabel("Median task R²")
            plt.title(f"{split_for_task_bars} median R² by task group")
            outputs["group_median_r2"] = savefig(
                f"{split_for_task_bars}_group_median_r2.png"
            )

            plt.figure(figsize=(12, 6))
            plt.bar(
                group_split_df["group"],
                group_split_df["median_pred_std_over_label_std"],
                color="tab:green",
            )
            plt.axhline(0.05, color="red", linestyle="--", linewidth=1)
            plt.axhline(0.25, color="orange", linestyle="--", linewidth=1)
            plt.axhline(1.0, color="black", linewidth=1)
            plt.xticks(rotation=45, ha="right")
            plt.ylabel("Median prediction std / label std")
            plt.title(
                f"{split_for_task_bars} median prediction variance ratio by group"
            )
            outputs["group_prediction_variance_ratio"] = savefig(
                f"{split_for_task_bars}_group_prediction_variance_ratio.png"
            )

    # -------------------------------------------------------------------------
    # 6. Task overlap heatmaps
    # -------------------------------------------------------------------------
    if overlap_counts_df is not None and overlap_counts_df.shape[0] > 0:
        matrix = overlap_counts_df.values.astype(float)

        plt.figure(figsize=heatmap_figsize)
        plt.imshow(matrix, aspect="auto")
        plt.colorbar(label="Number of compounds labelled for both tasks")
        plt.xticks(
            np.arange(overlap_counts_df.shape[1]),
            [wrap_label(c, width=16) for c in overlap_counts_df.columns],
            rotation=90,
            fontsize=6,
        )
        plt.yticks(
            np.arange(overlap_counts_df.shape[0]),
            [wrap_label(i, width=16) for i in overlap_counts_df.index],
            fontsize=6,
        )
        plt.title("Task overlap counts")
        outputs["task_overlap_counts_heatmap"] = savefig(
            "task_overlap_counts_heatmap.png"
        )

    if overlap_fraction_df is not None and overlap_fraction_df.shape[0] > 0:
        matrix = overlap_fraction_df.values.astype(float)

        plt.figure(figsize=heatmap_figsize)
        plt.imshow(matrix, aspect="auto", vmin=0.0, vmax=1.0)
        plt.colorbar(label="Fraction of row task labels also labelled for column task")
        plt.xticks(
            np.arange(overlap_fraction_df.shape[1]),
            [wrap_label(c, width=16) for c in overlap_fraction_df.columns],
            rotation=90,
            fontsize=6,
        )
        plt.yticks(
            np.arange(overlap_fraction_df.shape[0]),
            [wrap_label(i, width=16) for i in overlap_fraction_df.index],
            fontsize=6,
        )
        plt.title("Task overlap fraction by row task")
        outputs["task_overlap_fraction_heatmap"] = savefig(
            "task_overlap_fraction_heatmap.png"
        )

    # -------------------------------------------------------------------------
    # 7. Gradient cosine heatmap
    # -------------------------------------------------------------------------
    if grad_cosine_df is not None and len(grad_cosine_df) > 0:
        if "task" in grad_cosine_df.columns:
            grad_matrix_df = grad_cosine_df.set_index("task")
        else:
            grad_matrix_df = grad_cosine_df.copy()

        matrix = grad_matrix_df.values.astype(float)

        plt.figure(figsize=heatmap_figsize)
        plt.imshow(matrix, aspect="auto", vmin=-1.0, vmax=1.0, cmap="coolwarm")
        plt.colorbar(label="Gradient cosine similarity")
        plt.xticks(
            np.arange(grad_matrix_df.shape[1]),
            [wrap_label(c, width=16) for c in grad_matrix_df.columns],
            rotation=90,
            fontsize=6,
        )
        plt.yticks(
            np.arange(grad_matrix_df.shape[0]),
            [wrap_label(i, width=16) for i in grad_matrix_df.index],
            fontsize=6,
        )
        plt.title("Task gradient cosine similarity on shared parameters")
        outputs["task_gradient_cosine_heatmap"] = savefig(
            "task_gradient_cosine_heatmap.png"
        )

        # Per-task conflict score: mean cosine to all other tasks.
        conflict_rows = []

        for task in grad_matrix_df.index:
            values = grad_matrix_df.loc[task].astype(float)
            values = values.drop(labels=[task], errors="ignore")
            values = values.replace([np.inf, -np.inf], np.nan)

            conflict_rows.append({
                "task": task,
                "mean_cosine_to_other_tasks": values.mean(skipna=True),
                "fraction_negative_cosines": (values < 0).mean(skipna=True),
            })

        conflict_df = pd.DataFrame(conflict_rows)
        conflict_path = os.path.join(plot_dir, "task_gradient_conflict_summary.csv")
        conflict_df.to_csv(conflict_path, index=False)
        outputs["task_gradient_conflict_summary_csv"] = conflict_path
        logger.info(f"Saved gradient conflict summary: {conflict_path}")

        conflict_df = conflict_df.sort_values(
            "mean_cosine_to_other_tasks",
            ascending=True,
        )

        plt.figure(figsize=bar_figsize)
        plt.barh(
            [wrap_label(t) for t in conflict_df["task"]],
            conflict_df["mean_cosine_to_other_tasks"],
            color=[
                "tab:red" if v < 0 else "tab:blue"
                for v in conflict_df["mean_cosine_to_other_tasks"]
            ],
        )
        plt.axvline(0.0, color="black", linewidth=1)
        plt.xlabel("Mean gradient cosine to other tasks")
        plt.ylabel("Task")
        plt.title("Per-task mean gradient compatibility")
        outputs["task_gradient_mean_cosine"] = savefig(
            "task_gradient_mean_cosine_by_task.png"
        )

    # -------------------------------------------------------------------------
    # 8. Predicted vs observed scatter grids for worst and best validation tasks
    # -------------------------------------------------------------------------
    if metrics_df is not None:
        split_metrics = metrics_df[
            metrics_df["split"] == scatter_split
        ].copy()

        pred_df = read_csv_if_exists(f"{scatter_split}_predictions_scaled.csv")
        label_df = read_csv_if_exists(f"{scatter_split}_labels_scaled.csv")

        if (
            len(split_metrics) > 0 and
            pred_df is not None and
            label_df is not None
        ):
            split_metrics = split_metrics[
                split_metrics["n"].fillna(0) >= 2
            ].copy()

            worst_tasks = (
                split_metrics.sort_values("r2", ascending=True)
                .head(worst_n_scatter_tasks)["task"]
                .tolist()
            )

            best_tasks = (
                split_metrics.sort_values("r2", ascending=False)
                .head(best_n_scatter_tasks)["task"]
                .tolist()
            )

            def plot_scatter_grid(tasks, filename, title):
                if not tasks:
                    return None

                n_tasks = len(tasks)
                ncols = 3
                nrows = int(math.ceil(n_tasks / ncols))

                plt.figure(figsize=scatter_figsize)

                for i, task in enumerate(tasks):
                    ax = plt.subplot(nrows, ncols, i + 1)

                    pred_col = f"pred_{task}"
                    label_col = f"label_{task}"

                    if pred_col not in pred_df.columns or label_col not in label_df.columns:
                        ax.set_axis_off()
                        continue

                    merged = label_df[["row_index", label_col]].merge(
                        pred_df[["row_index", pred_col]],
                        on="row_index",
                        how="inner",
                    )

                    x = merged[label_col].values.astype(float)
                    y = merged[pred_col].values.astype(float)

                    mask = np.isfinite(x) & np.isfinite(y)
                    x = x[mask]
                    y = y[mask]

                    if len(x) == 0:
                        ax.set_axis_off()
                        continue

                    ax.scatter(x, y, alpha=0.5, s=12)

                    min_v = np.nanmin([np.min(x), np.min(y)])
                    max_v = np.nanmax([np.max(x), np.max(y)])

                    ax.plot([min_v, max_v], [min_v, max_v], color="black", linewidth=1)

                    task_metric = split_metrics[split_metrics["task"] == task]

                    if len(task_metric) > 0:
                        r2 = task_metric["r2"].iloc[0]
                        ratio = task_metric["pred_std_over_label_std"].iloc[0]
                        n = task_metric["n"].iloc[0]
                        subtitle = f"R²={r2:.2f}, ratio={ratio:.2f}, n={int(n)}"
                    else:
                        subtitle = ""

                    ax.set_title(
                        f"{wrap_label(task, width=24)}\n{subtitle}",
                        fontsize=8,
                    )
                    ax.set_xlabel("Scaled label")
                    ax.set_ylabel("Scaled prediction")

                plt.suptitle(title)
                return savefig(filename)

            outputs["worst_task_scatter_grid"] = plot_scatter_grid(
                worst_tasks,
                f"{scatter_split}_worst_tasks_pred_vs_label.png",
                f"{scatter_split} predicted vs observed, worst tasks by R²",
            )

            outputs["best_task_scatter_grid"] = plot_scatter_grid(
                best_tasks,
                f"{scatter_split}_best_tasks_pred_vs_label.png",
                f"{scatter_split} predicted vs observed, best tasks by R²",
            )

    # -------------------------------------------------------------------------
    # 9. Gradient batch counts
    # -------------------------------------------------------------------------
    if grad_counts_df is not None and len(grad_counts_df) > 0:
        counts_df = grad_counts_df.sort_values(
            "n_batches_with_labels",
            ascending=True,
        )

        plt.figure(figsize=bar_figsize)
        plt.barh(
            [wrap_label(t) for t in counts_df["task"]],
            counts_df["n_batches_with_labels"],
            color="tab:cyan",
        )
        plt.xlabel("Number of gradient batches with labels")
        plt.ylabel("Task")
        plt.title(
            "Task label availability during gradient cosine calculation\n"
            "Low values mean gradient estimates are less reliable"
        )
        outputs["gradient_batch_counts"] = savefig(
            "task_gradient_batch_counts.png"
        )

    manifest_path = os.path.join(plot_dir, "diagnostic_plot_outputs.json")

    with open(manifest_path, "w") as f:
        json.dump(outputs, f, indent=2)

    logger.info(f"Saved diagnostic plot manifest: {manifest_path}")
    logger.info("Finished diagnostic plotting")
    logger.info("=" * 80)

    context["diagnostic_plot_outputs"] = outputs

    return context