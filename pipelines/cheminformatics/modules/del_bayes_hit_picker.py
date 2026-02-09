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
import duckdb
import warnings

import polars as pl
import pandas as pd
import numpy as np
import pymc as pm
import arviz as az
from pathlib import Path
from itertools import combinations

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger
from modules.utils.plot_del_hits import plot_del_hits, plot_squares_hits, reduce_hits_for_plotting
from modules.utils.umap import add_umap_clusters

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

logger.debug("[CPU LIMIT] Affinity set to cores:", p.cpu_affinity())

@register_task(
    "del_stream_aggregate_counts",
    category="DEL",
    description="Stream DEL CSV chunks and aggregate synthon-level counts (DuckDB)."
)
def del_stream_aggregate_counts(config: dict, data: dict) -> dict:
    params = config.get("del_stream_aggregate_counts", {})
    group_cols = params.get("group_cols", ["NsynthonID", "lib", "condition"])
    count_col = params.get("count_col", "copy")
    physchem_cols = params.get("physchem_cols", ["RCount", "ARCount", "Fsp3", "logP", "MW"])
    min_pre = params.get("min_pre_count", 0)

    df = data.get("df")
    if df is None or len(df) == 0:
        return {"df": pd.DataFrame()}

    df[count_col] = pd.to_numeric(df[count_col], errors="coerce").fillna(0)
    for col in physchem_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # DuckDB in-memory table
    con = duckdb.connect(":memory:")
    con.register("del_df", df)

    agg_exprs = [
            "SUM(copy) AS post_count",
            "AVG(enrichment) AS enrichment_mean",  # for reporting only
            "COUNT(*) AS n_compounds",
        ]
    for col in physchem_cols:
        if col in df.columns:
            agg_exprs.append(f"AVG({col}) AS {col}_mean")

    select_str = ", ".join(group_cols + agg_exprs)
    group_str = ", ".join(group_cols)
    query = f"SELECT {select_str} FROM del_df GROUP BY {group_str}"
    out_df = con.execute(query).df()

    if min_pre > 0 and "total_count" in out_df.columns:
        out_df = out_df[out_df["total_count"] >= min_pre]

    # Propagate parquet path set by workflow runner
    parquet_path = data.get("output_path")  # workflow runner sets this
    if parquet_path:
        data["parquet_path"] = Path(parquet_path)

    return {"df": out_df, "parquet_path": data.get("parquet_path")}

# # Top-level function for ADVI chains to avoid multiprocessing pickling issues
# def run_advi_chain(seed, y_obs, cond_idx, syn_idx, cond_codes, syn_codes, overdisp, n_advi, draws, threads_per_chain):

#     # Limit threads per chain
#     for var in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
#                 "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"]:
#         os.environ[var] = str(threads_per_chain)

#     with pm.Model() as model:
#         alpha_cond = pm.Normal("alpha_cond", mu=0.0, sigma=1.0, shape=len(cond_codes))
#         sigma_syn = pm.HalfNormal("sigma_syn", sigma=1.0)
#         beta_syn = pm.Normal("beta_syn", mu=0.0, sigma=sigma_syn, shape=len(syn_codes))
#         log_mu = alpha_cond[cond_idx] + beta_syn[syn_idx]
#         mu = pm.math.exp(log_mu)
#         pm.NegativeBinomial("obs", mu=mu, alpha=overdisp, observed=y_obs)

#         approx = pm.fit(n=5000, 
#                         method="advi", 
#                         obj_optimizer=pm.adam(learning_rate=0.01),
#                         progressbar=True, 
#                         random_seed=int(seed)
#                         )
#         return approx.sample(draws=draws)
    
#     final_loss = approx.hist[-1]
#     initial_loss = approx.hist[0]

#     logger.info(
#         "[del_bayesian_model] ADVI loss: %.3e → %.3e",
#         initial_loss,
#         final_loss,
#     )

@register_task(
    "del_bayesian_model",
    category="DEL",
    description="Fit hierarchical Bayesian model to synthon enrichment per condition using ADVI with low-count filtering."
)
def del_bayesian_model(config: dict, data: dict) -> dict:

    import time
    import numpy as np
    import pandas as pd
    import pymc as pm
    import arviz as az
    import polars as pl
    import psutil
    from pathlib import Path

    params = config.get("del_bayesian_model", {})
    model_cfg = params.get("model", {})

    # ----------------------------
    # Resource check
    # ----------------------------
    avail_gb = psutil.virtual_memory().available / 1e9
    logger.info("[del_bayesian_model] Available RAM at start: %.1f GB", avail_gb)

    # ----------------------------
    # Load dataframe
    # ----------------------------
    df = data.get("df")
    if df is None or len(df) == 0:
        parquet_path = data.get("input_file")
        if parquet_path and Path(parquet_path).exists():
            df = pl.read_parquet(parquet_path)
        else:
            raise ValueError("[del_bayesian_model] Input dataframe empty and no fallback found.")

    # ----------------------------
    # Normalise schema
    # ----------------------------
    if "copy" not in df.columns:
        if "post_count" in df.columns:
            df = df.with_columns(pl.col("post_count").alias("copy")) if isinstance(df, pl.DataFrame) else df.assign(copy=df["post_count"])
        else:
            raise ValueError("[del_bayesian_model] No count column found (copy/post_count).")

    if "enrichment" not in df.columns:
        if "enrichment_mean" in df.columns:
            df = df.with_columns(pl.col("enrichment_mean").alias("enrichment")) if isinstance(df, pl.DataFrame) else df.assign(enrichment=df["enrichment_mean"])
        else:
            raise ValueError("[del_bayesian_model] No enrichment column found (enrichment/enrichment_mean).")

    if isinstance(df, pl.DataFrame):
        df = df.to_pandas()
    else:
        df = df.copy()

    # ----------------------------
    # Column validation
    # ----------------------------
    required = {"NsynthonID", "condition", "copy", "enrichment"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"[del_bayesian_model] Missing required columns: {missing}")

    # ----------------------------
    # Numeric coercion
    # ----------------------------
    df["copy"] = pd.to_numeric(df["copy"], errors="coerce").fillna(0).astype(int)
    df["enrichment"] = pd.to_numeric(df["enrichment"], errors="coerce")
    df = df.dropna(subset=["enrichment"])

    # ----------------------------
    # Apply min_pre_count (copy threshold)
    # ----------------------------
    min_pre = params.get("min_pre_count", 5)  # Increase for debugging stability
    df = df[df["copy"] >= min_pre]

    if df.empty:
        raise ValueError("[del_bayesian_model] No rows left after min_pre_count filtering.")

    # ----------------------------
    # Collapse to synthon × condition
    # ----------------------------
    model_df = (
        df.groupby(["NsynthonID", "condition"], as_index=False)
        .agg(
            enrichment_mean=("enrichment", "mean"),
            n_compounds=("copy", "count"),
        )
    )

    logger.info(
        "[del_bayesian_model] Collapsed %d rows → %d synthon-condition rows",
        len(df),
        len(model_df),
    )

    # ----------------------------
    # Response: log enrichment with clipping
    # ----------------------------
    y = model_df["enrichment_mean"].to_numpy(dtype=float)
    y = np.clip(y, 1e-6, 1000.0)  # avoid log(0) and huge spikes
    log_y = np.log(y)
    weights = np.sqrt(model_df["n_compounds"].to_numpy())  # weight by compound count

    # ----------------------------
    # Encode indices
    # ----------------------------
    syn_codes, syn_idx = np.unique(model_df["NsynthonID"], return_inverse=True)
    cond_codes, cond_idx = np.unique(model_df["condition"], return_inverse=True)

    # ----------------------------
    # Model parameters
    # ----------------------------
    n_advi = model_cfg.get("advi_iterations", 3000)
    draws = model_cfg.get("draws", 500)
    random_seed = model_cfg.get("random_seed", 42)

    logger.info("[del_bayesian_model] ADVI start: iterations=%d, synthons=%d, conditions=%d",
                n_advi, len(syn_codes), len(cond_codes))

    # ----------------------------
    # Bayesian model with partial pooling across synthons
    # ----------------------------
    start_time = time.time()
    with pm.Model() as model:

        n_cond = len(cond_codes)
        n_obs = len(model_df)
        
        # --- HIERARCHICAL PARTIAL POOLING ---
        # Number of unique synthons
        n_syn = int(syn_idx.max()) + 1

        # Synthon-level effects (partial pooling anchor)
        beta_syn = pm.Normal(
            "beta_syn",
            mu=0.0,
            sigma=1.0,
            shape=n_syn
        )

        # Per-condition shrinkage
        sigma_cond = pm.HalfNormal(
            "sigma_cond",
            sigma=1.0,
            shape=n_cond
        )

        # Synthon × condition effects (one per row in model_df)
        beta_syn_cond = pm.Normal(
            "beta_syn_cond",
            mu=beta_syn[syn_idx], 
            sigma=sigma_cond[cond_idx],
            shape=n_obs
        )

        # Observation noise
        sigma_obs = pm.HalfNormal("sigma_obs", sigma=1.0)

        # Likelihood
        pm.Normal(
            "obs",
            mu=beta_syn_cond,
            sigma=sigma_obs / weights,
            observed=log_y
        )

        approx = pm.fit(
            n=n_advi,
            method="advi",
            obj_optimizer=pm.adam(learning_rate=0.01),
            progressbar=True,
            random_seed=random_seed
        )

        trace = approx.sample(draws=draws)

    elapsed = (time.time() - start_time) / 60
    logger.info("[del_bayesian_model] ADVI finished in %.1f minutes", elapsed)

    # ----------------------------
    # Posterior summaries
    # ----------------------------
    summary_df = model_df[["NsynthonID", "condition"]].copy()

    beta_post = trace.posterior["beta_syn_cond"]  # (chain, draw, obs)

    # Mean
    summary_df["beta_mean"] = (
        beta_post.mean(dim=("chain", "draw"))
        .values
    )

    # HDI
    hdi = az.hdi(beta_post, hdi_prob=0.95)
    summary_df["beta_hdi_lower"] = hdi["beta_syn_cond"].sel(hdi="lower").values
    summary_df["beta_hdi_upper"] = hdi["beta_syn_cond"].sel(hdi="higher").values

    # Probability of activity
    summary_df["p_active"] = (
        (beta_post > 0)
        .mean(dim=("chain", "draw"))
        .values
    )

    summary_df["p_active_cond"] = summary_df["p_active"]

    logger.debug(
        "[del_bayesian_model] p_active quantiles: %s",
        np.quantile(summary_df["p_active"], [0.5, 0.8, 0.9, 0.95])
    )

    # ----------------------------
    # Merge physchem columns (mean per synthon-condition)
    # ----------------------------
    base_physchem_cols = [
        "MW","logP","PSA","ROTN","HBD","HBA",
        "RCount","ARCount","Fsp3","HeavyAcount"
    ]

    # Log current df columns for debugging
    logger.debug("[del_bayesian_model] Columns in df before physchem merge: %s", df.columns.tolist())

    physchem_cols = [c for c in df.columns if any(c == bc or c == f"{bc}_mean" for bc in base_physchem_cols)]

    logger.debug("[del_bayesian_model] physchem_cols present in df: %s", physchem_cols)

    if physchem_cols:
        physchem_df = df.groupby(["NsynthonID","condition"], as_index=False).agg(
            {c: "mean" for c in physchem_cols}
        )
        # HACK: rename merged columns with "_mean" to avoid conflicts
        physchem_df = physchem_df.rename(columns={c: f"{c}_mean" for c in physchem_cols})

        summary_df = summary_df.merge(
            physchem_df,
            on=["NsynthonID","condition"],
            how="left"
        )

        physchem_cols = [f"{c}_mean" for c in physchem_cols]
    else:
        physchem_cols = []

    logger.debug("[del_bayesian_model] Columns in summary_df after physchem merge: %s", summary_df.columns.tolist())
    logger.debug("[del_bayesian_model] Physchem columns returned: %s", physchem_cols)
    logger.debug("[del_bayesian_model] Non-NaN counts per physchem column:")
    for c in physchem_cols:
        col_mean = f"{c}_mean"
        if col_mean in summary_df.columns:
            logger.info("  %s: %d non-NaN", col_mean, summary_df[col_mean].notna().sum())

    return {"df": summary_df, "trace": trace, "physchem_cols": physchem_cols}

@register_task(
    "del_hit_picker",
    category="DEL",
    description="Select DEL synthons based on Bayesian posteriors and generate plots."
)
def del_hit_picker(config: dict, data: dict) -> dict:
    
    params = config.get("del_hit_picker", {})
    df = data.get("df")  # Model output from del_bayesian_model

    if df is None or len(df) == 0:
        raise ValueError("[del_hit_picker] Input dataframe is empty.")

    # Convert to pandas if needed
    if isinstance(df, pl.DataFrame):
        df = df.to_pandas()
    else:
        df = df.copy()

    # ----------------------------
    # Identify physchem columns
    # ----------------------------
    physchem_cols = data.get("physchem_cols", None)
    
    logger.debug("[del_hit_picker] physchem_cols from data dict: %s", physchem_cols)

    # Keep only columns that actually exist and are numeric
    if physchem_cols:
        physchem_cols = [c for c in physchem_cols if f"{c}_mean" in df.columns]
    
    logger.debug("[del_hit_picker] physchem_cols after existence check in df: %s", physchem_cols)

    # fallback
    if not physchem_cols:
        physchem_cols = [
            c for c in df.columns
            if c not in {"NsynthonID", "condition", "beta_mean", "beta_hdi_lower", "beta_hdi_upper", "p_active", "p_active_cond"}
            and pd.api.types.is_numeric_dtype(df[c])
            and df[c].notna().any()
        ]
    logger.debug("[del_hit_picker] Final physchem_cols used for plotting: %s", physchem_cols)
    logger.debug("[del_hit_picker] Non-NaN counts per physchem column:")

    for c in physchem_cols:
        logger.debug("  %s: %d non-NaN", c, df[c].notna().sum())

    logger.debug("[del_hit_picker] Physchem columns available for plotting: %s", physchem_cols)

    # ----------------------------
    # Apply Bayesian filters
    # ----------------------------

    posterior_cutoff = params.get("posterior_cutoff", 0.95)
    require_hdi_positive = params.get("require_hdi_positive", True)
    physchem_filters = params.get("physchem_filters", {})

    n_start = len(df)
    logger.info("[del_hit_picker] Starting with %d rows", n_start)

    mask = df["p_active"] >= posterior_cutoff
    logger.info(
        "[del_hit_picker] p_active >= %.3f: %d / %d pass",
        posterior_cutoff,
        mask.sum(),
        n_start,
    )

    if require_hdi_positive:
        mask_hdi = df["beta_hdi_lower"] > 0.0
        logger.info(
            "[del_hit_picker] beta_hdi_lower > 0: %d / %d pass",
            (mask & mask_hdi).sum(),
            n_start,
        )
        mask &= mask_hdi

    # ----------------------------
    # Apply physchem filters if provided
    # ----------------------------

    for col, (low, high) in physchem_filters.items():
        if col in df.columns:
            mask_phys = df[col].between(low, high)
            logger.info(
                "[del_hit_picker] %s in [%.3f, %.3f]: %d / %d pass",
                col,
                low,
                high,
                (mask & mask_phys).sum(),
                n_start
            )
            mask &= mask_phys

    # ----------------------------
    # Filtered hits
    # ----------------------------

    hits_df = df[mask].copy()

    if hits_df.empty:
        logger.warning("[del_hit_picker] No hits found after filtering.")
        top_hits_df = pd.DataFrame(columns=[
            "NsynthonID", "top_condition", "top_beta_mean",
            "top_beta_hdi_lower", "top_beta_hdi_upper", "top_p_active"
        ])
        return {"df": hits_df, "top_hits_df": top_hits_df, "plot_files": {}}

    # ----------------------------
    # Sort hits
    # ----------------------------

    sort_by = params.get("sort_by", ["p_active", "beta_mean"])
    ascending = params.get("ascending", [False, False])
    hits_df = hits_df.sort_values(by=sort_by, ascending=ascending)

    # ----------------------------
    # Top-hit-per-synthon summary
    # ----------------------------

    top_hits_df = (
        hits_df.sort_values(["p_active", "beta_mean"], ascending=[False, False])
        .groupby("NsynthonID", as_index=False)
        .first()
        .rename(columns={
            "condition": "top_condition",
            "beta_mean": "top_beta_mean",
            "beta_hdi_lower": "top_beta_hdi_lower",
            "beta_hdi_upper": "top_beta_hdi_upper",
            "p_active": "top_p_active"
        })
    )

    # ----------------------------
    # Save CSVs
    # ----------------------------

    out_dir = Path(params.get("output", {}).get("directory", "outputs/del_hits"))
    out_dir.mkdir(parents=True, exist_ok=True)
    full_csv_path = out_dir / params.get("output", {}).get("filename", "wuxi_del_hits.csv")
    top_csv_path = out_dir / params.get("output", {}).get("top_hits_filename", "wuxi_del_top_hits.csv")

    hits_df.to_csv(full_csv_path, index=False)
    top_hits_df.to_csv(top_csv_path, index=False)

    logger.info("[del_hit_picker] Full hits CSV saved: %s", str(full_csv_path))
    logger.info("[del_hit_picker] Top-hit summary CSV saved: %s", str(top_csv_path))

    logger.info(
        "[del_hit_picker] Running clustering step (UMAP)"
    )

    # ----------------------------
    # UMAP clustering for plotting
    # ----------------------------

    umap_features = [
        "beta_mean",
        "p_active",
    ]

    umap_features += [c for c in physchem_cols if c in hits_df.columns]

    hits_df = add_umap_clusters(
        hits_df,
        feature_cols=umap_features,
        n_neighbors=30,
        min_dist=0.15,
        cluster_method="hdbscan",
    )

    # ----------------------------
    # Reduce hits for plotting only
    # ----------------------------
    
    max_plot_hits = params.get("max_plot_hits", 1000) # TODO: ADD THIS TO YAML

    reduced_hits_df = reduce_hits_for_plotting(
        hits_df,
        max_hits=max_plot_hits
    )

    reduced_csv_path = out_dir / "wuxi_del_hits_reduced_umap.csv"
    reduced_hits_df.to_csv(reduced_csv_path, index=False)

    logger.info(
        "[del_hit_picker] Reduced (UMAP) hits CSV saved: %s",
        str(reduced_csv_path)
    )

    # ----------------------------
    # Generate plots
    # ----------------------------

    plot_dir = params.get("output", {}).get("directory", "outputs/plots")
    Path(plot_dir).mkdir(parents=True, exist_ok=True)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=FutureWarning)

        # Ensure physchem columns are numeric
        for c in physchem_cols:
            hits_df[c] = pd.to_numeric(hits_df[c], errors="coerce")

        plot_files = plot_del_hits(
            hits_df=reduced_hits_df,   # reduced df for visuals
            full_hits_df=hits_df,      # full df for UpSet
            output_dir=plot_dir,
            top_n=params.get("top_n_per_condition", 10),
            physchem_cols=physchem_cols,
            heatmap_top_n=50
        )

    plot_files_str = {k: str(v) for k, v in plot_files.items()}
    logger.debug("[del_hit_picker] Plots generated: %s", plot_files_str)

    return {
        "df": hits_df,
        "top_hits_df": top_hits_df,
        "reduced_hits_df": reduced_hits_df,
        "plot_files": plot_files_str
    }

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
    import pandas as pd
    from pathlib import Path

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