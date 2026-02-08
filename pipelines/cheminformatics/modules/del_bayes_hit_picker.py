# before importing pymc, somewhere in your startup code
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

# # 2. Limit BLAS / OpenMP / NumExpr threads (must happen before numpy/pymc import)
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

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger
from modules.utils.plot_del_hits import plot_del_hits

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
    # Normalize schema
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
    # Optional: Debug subset (top synthons only)
    # ----------------------------
    # debug_top = params.get("debug_top_synthons", 200)
    # if debug_top > 0:
    #     top_synthons = df.groupby("NsynthonID")["copy"].sum().nlargest(debug_top).index
    #     df = df[df["NsynthonID"].isin(top_synthons)]
    #     logger.info("[del_bayesian_model] Debug subset: %d synthons selected", len(top_synthons))

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
    # Bayesian model
    # ----------------------------
    start_time = time.time()
    with pm.Model() as model:
        alpha_cond = pm.Normal("alpha_cond", mu=0.0, sigma=1.0, shape=len(cond_codes))
        sigma_syn = pm.HalfNormal("sigma_syn", sigma=1.0)
        beta_syn = pm.Normal("beta_syn", mu=0.0, sigma=sigma_syn, shape=len(syn_codes))
        sigma_obs = pm.HalfNormal("sigma_obs", sigma=1.0)

        mu = alpha_cond[cond_idx] + beta_syn[syn_idx]

        pm.Normal("obs", mu=mu, sigma=sigma_obs / weights, observed=log_y)

        approx = pm.fit(n=n_advi, method="advi",
                        obj_optimizer=pm.adam(learning_rate=0.01),
                        progressbar=True, random_seed=random_seed)
        trace = approx.sample(draws=draws)

    elapsed = (time.time() - start_time) / 60
    logger.info("[del_bayesian_model] ADVI finished in %.1f minutes", elapsed)

    # ----------------------------
    # Posterior summaries (corrected per synthon × condition)
    # ----------------------------
    
    summary_df = model_df[["NsynthonID", "condition"]].copy()

    # Posterior arrays
    beta_post = trace.posterior["beta_syn"]      # shape: (chains, draws, n_syn)
    alpha_post = trace.posterior["alpha_cond"]   # shape: (chains, draws, n_cond)

    # Means
    beta_mean = beta_post.mean(dim=("chain","draw")).values  # synthon-level
    alpha_mean = alpha_post.mean(dim=("chain","draw")).values  # condition-level

    # HDIs
    beta_hdi = az.hdi(beta_post, hdi_prob=0.95)["beta_syn"].values
    alpha_hdi = az.hdi(alpha_post, hdi_prob=0.95)["alpha_cond"].values

    # Combine for synthon × condition
    summary_df["beta_mean"] = beta_mean[syn_idx] + alpha_mean[cond_idx]
    summary_df["beta_hdi_lower"] = beta_hdi[syn_idx,0] + alpha_hdi[cond_idx,0]
    summary_df["beta_hdi_upper"] = beta_hdi[syn_idx,1] + alpha_hdi[cond_idx,1]

    # Synthon-level p_active (correct)
    p_active = (beta_post > 0).mean(dim=("chain","draw")).values
    summary_df["p_active"] = p_active[syn_idx]


    # ----------------------------
    # Merge physchem (mean per synthon-condition)
    # ----------------------------
    physchem_cols = [
        "MW","logP","PSA","ROTN","HBD","HBA",
        "RCount","ARCount","Fsp3","HeavyAcount"
    ]
    physchem_cols = [c for c in physchem_cols if c in df.columns]

    if physchem_cols:
        physchem_df = df.groupby(["NsynthonID","condition"], as_index=False).agg({c:"mean" for c in physchem_cols})
        summary_df = summary_df.merge(physchem_df, on=["NsynthonID","condition"], how="left")

    return {"df": summary_df, "trace": trace}




@register_task(
    "del_hit_picker",
    category="DEL",
    description="Select DEL synthons based on Bayesian posteriors and generate plots."
)
def del_hit_picker(config: dict, data: dict) -> dict:
    
    params = config.get("del_hit_picker", {})
    df = data.get("df")  # Model output
    physchem_df = data.get("physchem_df")  # aggregated physchem data

    if df is None or len(df) == 0:
        raise ValueError("[del_hit_picker] Input dataframe is empty.")

    # Convert to pandas if needed
    if isinstance(df, pl.DataFrame):
        df = df.to_pandas()

    if physchem_df is not None and isinstance(physchem_df, pl.DataFrame):
        physchem_df = physchem_df.to_pandas()

    # Merge physchem features onto posterior table by synthon and condition
    physchem_cols = []
    if physchem_df is not None:
        # Keep only numeric physchem columns, exclude counts
        physchem_cols = [
            c for c in physchem_df.columns
            if c not in {"NsynthonID", "condition", "total_count", "n_compounds",
                         "beta_mean", "beta_hdi_lower", "beta_hdi_upper", "p_active"} 
            and pd.api.types.is_numeric_dtype(physchem_df[c])
        ]

        if physchem_cols:
            merge_cols = ["NsynthonID", "condition"] + physchem_cols
            df = df.merge(physchem_df[merge_cols], on=["NsynthonID", "condition"], how="left")

            # NaN check
            physchem_cols = [c for c in physchem_cols if c in df.columns and df[c].notna().any()]

    # ----------------------------
    # Apply Bayesian / physchem filters 
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

    for col, (low, high) in physchem_filters.items():
        if col in df.columns:
            mask_phys = df[col].between(low, high)
            logger.info(
                "[del_hit_picker] %s in [%.3f, %.3f]: %d / %d pass",
                col,
                low,
                high,
                (mask & mask_phys).sum(),
                n_start,
            )
            mask &= mask_phys

    hits_df = df[mask].copy()  # <-- filtered df, already has per-condition p_active

    # Check if no hits, gracefully fail and warn
    if hits_df.empty:
        logger.warning("[del_hit_picker] No hits found after filtering.")
        
        top_hits_df = pd.DataFrame(
            columns=[
                "NsynthonID", "top_condition", "top_beta_mean",
                "top_beta_hdi_lower", "top_beta_hdi_upper", "top_p_active"
            ]
        )
        
        return {
            "df": hits_df,          # empty
            "top_hits_df": top_hits_df,  # empty
            "plot_files": {}        # no plots
        }

    # Sorting 
    sort_by = params.get("sort_by", ["p_active", "beta_mean"])
    ascending = params.get("ascending", [False, False])
    hits_df = hits_df.sort_values(by=sort_by, ascending=ascending)

    # Top-hit-per-synthon summary 
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

    # Save CSVs 
    out_dir = Path(params.get("output", {}).get("directory", "outputs/del_hits"))
    out_dir.mkdir(parents=True, exist_ok=True)
    full_csv_path = out_dir / params.get("output", {}).get("filename", "wuxi_del_hits.csv")
    top_csv_path = out_dir / params.get("output", {}).get("top_hits_filename", "wuxi_del_top_hits.csv")

    # --- NO recomputation of p_active needed ---
    # hits_df["p_active"] already contains per-condition values used in plots

    hits_df.to_csv(full_csv_path, index=False)
    top_hits_df.to_csv(top_csv_path, index=False)

    logger.info("[del_hit_picker] Full hits CSV saved: %s", str(full_csv_path))
    logger.info("[del_hit_picker] Top-hit summary CSV saved: %s", str(top_csv_path))

    # Generate plots 
    plot_dir = params.get("output", {}).get("directory", "outputs/plots")
    Path(plot_dir).mkdir(parents=True, exist_ok=True)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=FutureWarning)

        # Ensure physchem columns are numeric
        for c in ["MW_mean", "logP_mean", "PSA_mean"]:
            if c in hits_df.columns:
                hits_df[c] = pd.to_numeric(hits_df[c], errors="coerce")

        valid_physchem_cols = [
            c for c in ["MW_mean", "logP_mean", "PSA_mean"]
            if c in hits_df.columns and hits_df[c].notna().any()
        ]

        logger.debug("Plotting with physchem columns: %s", valid_physchem_cols)
        logger.debug("Available columns: %s", list(hits_df.columns))

        plot_files = plot_del_hits(
            hits_df,
            output_dir=plot_dir,
            top_n=params.get("top_n_per_condition", 10),
            physchem_cols=valid_physchem_cols,
            heatmap_top_n = 50
        )

    plot_files_str = {k: str(v) for k, v in plot_files.items()}
    logger.info("[del_hit_picker] Plots generated: %s", plot_files_str)

    return {
        "df": hits_df,
        "top_hits_df": top_hits_df,
        "plot_files": plot_files_str
    }
