import random
from collections import defaultdict
from pathlib import Path
import joblib
import json

import numpy as np
import pandas as pd
import torch
import time

from tqdm import tqdm

# from itertools import combinations
# import matplotlib.cm as cm

import esm
from catboost import CatBoostRegressor

# from sklearn.model_selection import train_test_split
# from sklearn.decomposition import PCA
# from sklearn.neighbors import NearestNeighbors

# import umap
# import matplotlib.pyplot as plt
# import seaborn as sns
# from matplotlib import cm

from pipeline.task_registry import register_task
# from modules.utils.plotting import (
#     compute_multi_condition_pareto,
#     plot_round_radar,
#     generate_summary_plots,
#     plot_learning_curves_per_property,
# )
# from modules.utils.settings_logs import log_pipeline_config
# from modules.utils.eval_regression import save_evaluation, evaluate_with_bootstrap
# from modules.utils.eval_roundwise import run_roundwise_evaluation

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def make_unique(col, existing_cols):
    if col not in existing_cols:
        return col
    i = 1
    while f"{col}_{i}" in existing_cols:
        i += 1
    return f"{col}_{i}"

@register_task(
    "load_inference_dataset",
    category="VHH inference",
    description="Load sequences for property prediction"
)
def load_inference_dataset(config, context):

    df = pd.read_csv(config["csv_path"])

    seq_col = config.get("sequence_col", "sequence")
    expected_len = config.get("expected_len", 125)

    df = df.dropna(subset=[seq_col]).copy()

    df[seq_col] = df[seq_col].astype(str)
    df = df[
        df[seq_col].str.len() == expected_len
    ].reset_index(drop=True)

    context["df_inference"] = df
    context["sequence_col"] = seq_col

    logger.info(
        f"Loaded {len(df)} inference sequences"
    )

    return context

@register_task(
    "load_trained_model_bundle",
    category="VHH inference",
    description="Load PCA and CatBoost ensembles"
)
def load_trained_model_bundle(config, context):

    model_dir = Path(config["model_dir"])

    pca = joblib.load(
        model_dir / "pca.joblib"
    )

    with open(model_dir / "metadata.json") as f:
        meta = json.load(f)

    properties = meta["properties"]

    models = {}

    for prop in properties:

        ensemble = []

        i = 0
        while True:

            path = (
                model_dir /
                f"{prop.replace('/','_')}_model_{i}.cbm"
            )

            if not path.exists():
                break

            m = CatBoostRegressor()
            m.load_model(str(path))
            ensemble.append(m)

            i += 1

        models[prop] = ensemble

    context["pca"] = pca
    context["models"] = models
    context["properties"] = properties

    logger.info(
        f"Loaded models for {len(properties)} properties"
    )

    return context

@register_task(
    "run_model_inference",
    category="VHH inference",
    description="Predict properties on sequences"
)
def run_model_inference(config, context):

    df = context["df_inference"].copy()

    esm_model = context["esm_model"]
    batch_converter = context["batch_converter"]
    device = context["device"]

    pca = context["pca"]
    models = context["models"]
    properties = context["properties"]

    batch_size = config.get(
        "batch_size",
        32
    )

    seq_col = context["sequence_col"]

    @torch.no_grad()
    def embed_sequences(sequences):

        embs = []

        for i in tqdm(
            range(0,len(sequences),batch_size)
        ):

            batch = sequences[i:i+batch_size]

            data = [
                ("seq", s)
                for s in batch
            ]

            _,_,tokens = batch_converter(data)

            tokens = tokens.to(device)

            with torch.cuda.amp.autocast():
                out = esm_model(
                    tokens,
                    repr_layers=[33],
                    return_contacts=False
                )

            reps = out["representations"][33]

            for j,seq in enumerate(batch):

                emb = (
                    reps[j,1:len(seq)+1]
                    .mean(0)
                    .cpu()
                    .numpy()
                )

                embs.append(emb)

            del tokens,out,reps
            torch.cuda.empty_cache()

        return np.vstack(embs)


    logger.info("Embedding sequences...")
    X = embed_sequences(
        df[seq_col].tolist()
    )

    X_pca = pca.transform(X)

    logger.info("Running predictions...")

    for prop, ensemble in models.items():

        preds = np.column_stack(
            [
                m.predict(X_pca)
                for m in ensemble
            ]
        )

        mean_pred = preds.mean(axis=1)
        std_pred  = preds.std(axis=1)

        prefix = config.get("prediction_prefix", "PRED_")
        base_name = f"{prefix}{prop}"
        safe_name = make_unique(base_name, df.columns)

        df[safe_name] = mean_pred

        if config.get(
            "add_uncertainty",
            True
        ):
            df[f"{safe_name}_pred_std"] = std_pred
            df[f"{safe_name}_pred_lci"] = (
                mean_pred - 1.96*std_pred
            )
            df[f"{safe_name}_pred_uci"] = (
                mean_pred + 1.96*std_pred
            )

    out_csv = Path(config["output_csv"])

    out_csv.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(
        out_csv,
        index=False
    )

    logger.info(
        f"Saved predictions to {out_csv}"
    )

    return context