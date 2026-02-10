import os
import psutil

MAX_THREADS = min(16, psutil.cpu_count(logical=True))

for var in [
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
]:
    os.environ[var] = str(MAX_THREADS)

# ----------------------------
# GLOBAL CPU / THREAD LIMITS
# ----------------------------
MAX_CPUS = 128
# THREADS_PER_CHAIN = 8

# 1. Hard OS limit for the whole process
p = psutil.Process()
p.cpu_affinity(list(range(MAX_CPUS)))

# # 2. Limit BLAS / OpenMP / NumExpr threads
# os.environ["OMP_NUM_THREADS"] = str(THREADS_PER_CHAIN)
# os.environ["OPENBLAS_NUM_THREADS"] = str(THREADS_PER_CHAIN)
# os.environ["MKL_NUM_THREADS"] = str(THREADS_PER_CHAIN)
# os.environ["VECLIB_MAXIMUM_THREADS"] = str(THREADS_PER_CHAIN)
# os.environ["NUMEXPR_NUM_THREADS"] = str(THREADS_PER_CHAIN)

# ----------------------------
# IMPORTS AFTER THREAD LIMITS
# ----------------------------
# import duckdb
# import warnings

import polars as pl
import pandas as pd
# import numpy as np
# import pymc as pm
# import arviz as az
from pathlib import Path
from itertools import combinations

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger
from modules.utils.plot_squares_hits import plot_squares_hits
# from modules.utils.umap import add_umap_clusters

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

logger.debug("[CPU LIMIT] Affinity set to cores:", p.cpu_affinity())

@register_task(
    "del_squares_analysis",
    category="DEL",
    description="Perform proper squares analysis (SAR consistency) for DEL synthons."
)
def del_squares_analysis(config: dict, data: dict) -> dict:

    # ----------------------------
    # Config
    # ----------------------------

    params = config.get("del_squares_analysis", {})
    enrichment_col = params.get("enrichment_col", "enrichment_mean")  # column to use
    threshold = params.get("active_threshold", 1.5)  # enrichment threshold to call active
    out_dir = Path(params.get("output", {}).get("directory", "outputs/del_squares"))
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / params.get("output", {}).get("filename", "squares_counts.parquet")

    # ----------------------------
    # Get input data
    # ----------------------------

    df = data.get("df")
    if df is None or df.empty:
        parquet_path = data.get("input_file")
        if parquet_path and Path(parquet_path).exists():
            df = pl.read_parquet(parquet_path)
        else:
            raise ValueError("[del_squares_model] No input dataframe found and no valid input_file path.")
    
    # ----------------------------
    # Normalise schema
    # ----------------------------
    if "copy" not in df.columns:
        if "post_count" in df.columns:
            df = df.with_columns(pl.col("post_count").alias("copy")) if isinstance(df, pl.DataFrame) else df.assign(copy=df["post_count"])
        else:
            raise ValueError("[del_squares_model] No count column found (copy/post_count).")

    if "enrichment" not in df.columns:
        if "enrichment_mean" in df.columns:
            df = df.with_columns(pl.col("enrichment_mean").alias("enrichment")) if isinstance(df, pl.DataFrame) else df.assign(enrichment=df["enrichment_mean"])
        else:
            raise ValueError("[del_squares_model] No enrichment column found (enrichment/enrichment_mean).")

    if isinstance(df, pl.DataFrame):
        df = df.to_pandas()
    else:
        df = df.copy()

    # ----------------------------
    # Call active/inactive per condition
    # ----------------------------

    df["active"] = df[enrichment_col] >= threshold

    # ----------------------------
    # Enumerate squares per condition pair
    # ----------------------------

    conditions = df["condition"].unique()
    squares_list = []

    for c1, c2 in combinations(conditions, 2):
        df_c1 = df[df["condition"] == c1][["NsynthonID", "active"]].set_index("NsynthonID")
        df_c2 = df[df["condition"] == c2][["NsynthonID", "active"]].set_index("NsynthonID")

        # Synthon IDs active in both conditions
        active_c1 = set(df_c1[df_c1["active"]].index)
        active_c2 = set(df_c2[df_c2["active"]].index)
        square_synthon_ids = active_c1 & active_c2

        for syn in square_synthon_ids:
            squares_list.append(syn)

    # ----------------------------
    # Count nsquares per synthon
    # ----------------------------

    squares_count = pd.Series(squares_list).value_counts().reset_index()
    squares_count.columns = ["NsynthonID", "nsquares"]

    # Merge back per-condition enrichment info
    squares_df = df[['NsynthonID','condition',enrichment_col]].merge(
        squares_count, on='NsynthonID', how='right'
    )
    squares_df = squares_df.rename(columns={enrichment_col: 'enrichment_mean', 'nsquares': 'squares_score'})

    # ----------------------------
    # Save output
    # ----------------------------

    squares_df.to_parquet(out_file, index=False)

    return {
        "df": squares_df,
        "output_file": out_file
    }

@register_task(
    "del_squares_hit_picker",
    category="DEL",
    description="Pick hits from squares scores using thresholding and physchem filters, preserving per-condition info."
)
def del_squares_hit_picker(config: dict, data: dict) -> dict:

    params = config.get("del_squares_hit_picker", {})
    df = data.get("df")
    if df is None or df.empty:
        raise ValueError("[del_squares_hit_picker] No squares scores found.")

    score_cutoff = params.get("score_cutoff", 2.0)
    physchem_filters = params.get("physchem_filters", {})

    # ----------------------------
    # Filter by synthon-level score
    # ----------------------------

    hits_df = df[df["squares_score"] >= score_cutoff].copy()

    # ----------------------------
    # Apply physchem filters if provided
    # ----------------------------

    for col, (low, high) in physchem_filters.items():
        if col in hits_df.columns:
            hits_df = hits_df[hits_df[col].between(low, high)]

    # ----------------------------
    # Top-hit-per-synthon summary
    # ----------------------------

    top_hits_df = (
        hits_df.sort_values(["squares_score", "enrichment_mean"], ascending=[False, False])
        .groupby("NsynthonID", as_index=False)
        .first()
        .rename(columns={
            "condition": "top_condition",
            "enrichment_mean": "top_enrichment",
            "squares_score": "top_squares_score"
        })
    )

    # ----------------------------
    # Save CSVs
    # ----------------------------
    
    out_dir = Path(params.get("output", {}).get("directory", "outputs/del_squares_hits"))
    out_dir.mkdir(parents=True, exist_ok=True)
    hits_file = out_dir / params.get("output", {}).get("filename", "squares_hits.csv")
    top_file = out_dir / params.get("output", {}).get("top_hits_filename", "squares_top_hits.csv")

    hits_df.to_csv(hits_file, index=False)
    top_hits_df.to_csv(top_file, index=False)

    plot_files = plot_squares_hits(
        hits_df,
        top_n=50,
        enrichment_col="enrichment_mean",
        enrichment_threshold=1.5,
        thresholds=[1.5,2.0,2.5,3.0],
        output_dir=out_dir
    )

    return {
        "df": hits_df,
        "top_hits_df": top_hits_df,
        "output_files": {"full": hits_file, "top": top_file},
        "plot_files": {k: str(v) for k, v in plot_files.items()}
    }