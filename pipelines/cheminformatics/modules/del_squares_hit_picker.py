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

# 1. Hard OS limit for the whole process
p = psutil.Process()
p.cpu_affinity(list(range(MAX_CPUS)))
import umap

import polars as pl
import pandas as pd
from pathlib import Path
from itertools import combinations

from sklearn.preprocessing import StandardScaler

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger
from modules.utils.plot_squares_hits import plot_squares_hits

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

logger.debug("[CPU LIMIT] Affinity set to cores:", p.cpu_affinity())

@register_task(
    "del_squares_analysis",
    category="DEL",
    description="Perform Squares analysis for DEL synthons with optional UMAP clustering."
)
def del_squares_analysis(config: dict, data: dict) -> dict:

    # ----------------------------
    # Configs
    # ----------------------------

    params = config.get("del_squares_model", {})
    enrichment_col = params.get("enrichment_col", "enrichment_mean")
    threshold = params.get("active_threshold", 1.5)
    out_dir = Path(params.get("output", {}).get("directory", "outputs/del_squares"))
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / params.get("output", {}).get("filename", "squares_counts.parquet")

    # UMAP params
    use_umap = params.get("use_umap", True)
    umap_cols = params.get("umap_columns", ["squares_score", "enrichment_mean"])
    n_neighbors = params.get("umap_n_neighbors", 15)
    min_dist = params.get("umap_min_dist", 0.1)
    n_components = params.get("umap_n_components", 2)
    random_state = params.get("umap_random_state", 42)

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
        active_c1 = set(df_c1[df_c1["active"]].index)
        active_c2 = set(df_c2[df_c2["active"]].index)
        square_synthon_ids = active_c1 & active_c2
        squares_list.extend(square_synthon_ids)

    # ----------------------------
    # Count nsquares per synthon
    # ----------------------------
    squares_count = pd.Series(squares_list).value_counts().reset_index()
    squares_count.columns = ["NsynthonID", "squares_score"]

    squares_df = df[['NsynthonID','condition',enrichment_col]].drop_duplicates(subset=["NsynthonID", "condition"])
    squares_df = squares_df.merge(squares_count, on='NsynthonID', how='right')
    squares_df = squares_df.rename(columns={enrichment_col: 'enrichment_mean'})

    # Extract library from NsynthonID
    squares_df["library"] = squares_df["NsynthonID"].str.split("-", n=1).str[0]

    # Compute per-library weights
    lib_counts = squares_df["library"].value_counts()
    squares_df["umap_weight"] = 1.0 / lib_counts[squares_df["library"]].values

    # Sample hits for UMAP embedding
    n_umap_samples = min(len(squares_df), params.get("max_plot_hits", 5000))
    umap_input_df = squares_df.sample(
        n=n_umap_samples,
        weights="umap_weight",
        random_state=42
    ).reset_index(drop=True)

    # ----------------------------
    # Optional UMAP embedding
    # ----------------------------
    if use_umap and not squares_df.empty:
        umap_features = squares_df[umap_cols].fillna(0).values
        scaler = StandardScaler()
        umap_scaled = scaler.fit_transform(umap_features)
        reducer = umap.UMAP(
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            n_components=n_components,
            # random_state=random_state, # can't use with n_jobs
            init='random',
            verbose=True,
            n_jobs=20
        )
        embedding = reducer.fit_transform(umap_scaled)
        for i in range(n_components):
            squares_df[f"umap_{i+1}"] = embedding[:, i]
        squares_df["umap_cluster"] = embedding[:, 0].round().astype(int)

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
    out_dir = Path(params.get("output", {}).get("directory", "outputs/del_squares_hits"))
    out_dir.mkdir(parents=True, exist_ok=True)
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
    # Create reduced representative dataframe for plotting
    # ----------------------------
    if {"umap_1", "umap_2"}.issubset(hits_df.columns):

        # Extract library
        hits_df["library"] = hits_df["NsynthonID"].str.split("-", n=1).str[0]

        # Compute weights
        lib_counts = hits_df["library"].value_counts()
        hits_df["umap_weight"] = 1.0 / lib_counts[hits_df["library"]].values

        # Sample proportionally
        max_plot_hits = params.get("max_plot_hits", 2000)
        n_samples = min(len(hits_df), max_plot_hits)

        reduced_hits_df = (
            hits_df
            .sample(n=n_samples, weights="umap_weight", random_state=42)
            .reset_index(drop=True)
        )

        reduced_file = out_dir / "squares_reduced_UMAP.csv"
        reduced_hits_df.to_csv(reduced_file, index=False)

    else:
        logger.warning("UMAP columns not found; using full hits_df for plotting.")
        reduced_hits_df = hits_df
        reduced_file = None


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
    
    hits_file = out_dir / params.get("output", {}).get("filename", "squares_hits.csv")
    top_file = out_dir / params.get("output", {}).get("top_hits_filename", "squares_top_hits.csv")

    hits_df.to_csv(hits_file, index=False)
    top_hits_df.to_csv(top_file, index=False)

    plot_files = plot_squares_hits(
        df=reduced_hits_df,              # Reduced for plots
        full_df=hits_df,                 # Full df for other plots
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