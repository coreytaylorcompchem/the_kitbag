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

# Top-level function for ADVI chains to avoid multiprocessing pickling issues
def run_advi_chain(seed, y_obs, cond_idx, syn_idx, cond_codes, syn_codes, overdisp, n_advi, draws, threads_per_chain):

    # Limit threads per chain
    for var in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"]:
        os.environ[var] = str(threads_per_chain)

    with pm.Model() as model:
        alpha_cond = pm.Normal("alpha_cond", mu=0.0, sigma=1.0, shape=len(cond_codes))
        sigma_syn = pm.HalfNormal("sigma_syn", sigma=1.0)
        beta_syn = pm.Normal("beta_syn", mu=0.0, sigma=sigma_syn, shape=len(syn_codes))
        log_mu = alpha_cond[cond_idx] + beta_syn[syn_idx]
        mu = pm.math.exp(log_mu)
        pm.NegativeBinomial("obs", mu=mu, alpha=overdisp, observed=y_obs)

        approx = pm.fit(n=5000, 
                        method="advi", 
                        obj_optimizer=pm.adam(learning_rate=0.01),
                        progressbar=True, 
                        random_seed=int(seed)
                        )
        return approx.sample(draws=draws)
    
    final_loss = approx.hist[-1]
    initial_loss = approx.hist[0]

    logger.info(
        "[del_bayesian_model] ADVI loss: %.3e → %.3e",
        initial_loss,
        final_loss,
    )

@register_task(
    "del_bayesian_model",
    category="DEL",
    description="Fit hierarchical Bayesian model to synthon counts per condition using ADVI (optimized, multi-chain)."
)
def del_bayesian_model(config: dict, data: dict) -> dict:

    params = config.get("del_bayesian_model", {})
    model_cfg = params.get("model", {})

    avail_gb = psutil.virtual_memory().available / 1e9
    logger.info(
        "[del_bayesian_model] Available RAM at start: %.1f GB",
        avail_gb,
    )

    if avail_gb < 50:
        logger.warning(
            "[del_bayesian_model] Low available RAM — consider stopping early."
        )

    # Load dataframe (prefer Polars)
    df = data.get("df")
    if df is None or len(df) == 0:
        parquet_path = data.get("input_file")
        if parquet_path and Path(parquet_path).exists():
            df = pl.read_parquet(parquet_path)
        else:
            raise ValueError("[del_bayesian_model] Input dataframe empty and no fallback found.")

    if "post_count" not in df.columns:
        raise ValueError("[del_bayesian_model] Expected column 'post_count' in input dataframe.")

    # Convert to numeric
    if isinstance(df, pl.DataFrame):
        df = df.with_columns([pl.col("post_count").cast(pl.Int64)])
        df_pd = df.to_pandas()
    else:
        df["post_count"] = pd.to_numeric(df["post_count"], errors="coerce").fillna(0).astype(int)
        df_pd = df

    # ---- Collapse to synthon × condition (statistically equivalent, MUCH smaller)
    model_df = (
        df_pd
        .groupby(["NsynthonID", "condition"], as_index=False)
        .agg(post_count=("post_count", "sum"))
    )

    logger.info(
        "[del_bayesian_model] Collapsed %d rows → %d synthon-condition rows",
        len(df_pd),
        len(model_df),
    )

    logger.info(
        "[del_bayesian_model] Unique synthons=%d, conditions=%d, rows=%d",
        model_df["NsynthonID"].nunique(),
        model_df["condition"].nunique(),
        len(model_df),
    )

    # Encode synthons and conditions
    syn_codes, syn_idx = np.unique(model_df["NsynthonID"], return_inverse=True)
    cond_codes, cond_idx = np.unique(model_df["condition"], return_inverse=True)
    y_obs = model_df["post_count"].to_numpy()

    # Spike test
    spike_cfg = params.get("spike_test", {})
    if spike_cfg.get("enable", False):
        rng = np.random.default_rng(model_cfg.get("random_seed", 42))
        spike_idx = rng.choice(len(y_obs), size=min(spike_cfg.get("synthons", 10), len(y_obs)), replace=False)
        y_obs[spike_idx] *= spike_cfg.get("multiplier", 10)

    overdisp = model_cfg.get("overdispersion", 1.5)
    draws = model_cfg.get("draws", 500)
    chains = model_cfg.get("chains", 4)
    n_advi = model_cfg.get("advi_iterations", 10000) 
    random_seed = model_cfg.get("random_seed", 42)

    # Dynamic threads allocation
    # max_cpus = min(params.get("resources", {}).get("max_cpus", psutil.cpu_count(logical=True)), 64)
    # threads_per_chain = max(1, max_cpus // chains)
    # os.environ.update({
    #     "OMP_NUM_THREADS": str(threads_per_chain),
    #     "OPENBLAS_NUM_THREADS": str(threads_per_chain),
    #     "MKL_NUM_THREADS": str(threads_per_chain),
    #     "VECLIB_MAXIMUM_THREADS": str(threads_per_chain),
    #     "NUMEXPR_NUM_THREADS": str(threads_per_chain),
    # })

    # logger.info(f"[del_bayesian_model] ADVI with chains={chains}, threads_per_chain={threads_per_chain}, max_cpus={max_cpus}")

    # from multiprocessing import get_context

    # # Prepare args for each chain
    # args_list = [
    #     (seed, y_obs, cond_idx, syn_idx, cond_codes, syn_codes, overdisp, n_advi, draws, threads_per_chain)
    #     for seed in np.arange(chains) + random_seed
    # ]

    # # Use spawn context
    # with get_context("spawn").Pool(processes=chains) as pool:
    #     trace_list = pool.starmap(run_advi_chain, args_list)

    # # Combine traces
    # trace = az.concat(trace_list, dim="chain")

    zero_frac = (y_obs == 0).mean()
    logger.info(
        "[del_bayesian_model] Zero-count fraction: %.3f",
        zero_frac,
    )

    logger.info(
        "[del_bayesian_model] post_count quantiles:\n%s",
        pd.Series(y_obs).quantile([0, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0])
    )

    logger.info(
        "[del_bayesian_model] Running single ADVI fit."
    )

    import time

    with pm.Model() as model:
        alpha_cond = pm.Normal("alpha_cond", mu=0.0, sigma=1.0, shape=len(cond_codes))
        sigma_syn = pm.HalfNormal("sigma_syn", sigma=1.0)
        beta_syn = pm.Normal("beta_syn", mu=0.0, sigma=sigma_syn, shape=len(syn_codes))

        log_mu = alpha_cond[cond_idx] + beta_syn[syn_idx]
        mu = pm.math.exp(log_mu)

        pm.NegativeBinomial(
            "obs",
            mu=mu,
            alpha=overdisp,
            observed=y_obs,
        )

        start_time = time.time()

        logger.info(
            "[del_bayesian_model] ADVI start: iterations=%d, synthons=%d, conditions=%d",
            n_advi,
            len(syn_codes),
            len(cond_codes),
        )

        approx = pm.fit(
            n=n_advi,
            method="advi",
            obj_optimizer=pm.adam(learning_rate=0.01),
            progressbar=True,
            random_seed=random_seed,
        )

        elapsed = time.time() - start_time
        logger.info(
            "[del_bayesian_model] ADVI finished in %.1f minutes",
            elapsed / 60,
        )

        trace = approx.sample(draws=draws)

    # ----- Posterior statistics -----
    beta_post = trace.posterior["beta_syn"]
    p_active = (beta_post > 0).mean(dim=("chain", "draw")).values
    hdi = az.hdi(beta_post, hdi_prob=0.95)["beta_syn"].values
    beta_hdi_lower = hdi[..., 0] if hdi.ndim == 2 else hdi[:, 0]
    beta_hdi_upper = hdi[..., 1] if hdi.ndim == 2 else hdi[:, 1]
    beta_mean = beta_post.mean(dim=("chain", "draw")).values

    logger.info(
        "[del_bayesian_model] beta_mean summary:\n%s",
        pd.Series(beta_mean).describe()
    )

    summary_df = pd.DataFrame({
        "NsynthonID": df_pd["NsynthonID"],
        "condition": df_pd["condition"],
        "beta_mean": beta_mean[syn_idx],
        "beta_hdi_lower": beta_hdi_lower[syn_idx],
        "beta_hdi_upper": beta_hdi_upper[syn_idx],
        "p_active": p_active[syn_idx],
    })

    # Merge physchem features
    physchem_cols = [
        c for c in df_pd.columns
        if c.endswith("_mean") or c in ["MW","logP","PSA","ROTN","HBD","HBA","RCount","ARCount","Fsp3","HeavyAcount"]
    ]
    if physchem_cols:
        summary_df = summary_df.merge(
            df_pd[["NsynthonID", "condition"] + physchem_cols],
            on=["NsynthonID", "condition"],
            how="left"
        )

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

    # Apply Bayesian / physchem filters 
    posterior_cutoff = params.get("posterior_cutoff", 0.95)
    require_hdi_positive = params.get("require_hdi_positive", True)
    physchem_filters = params.get("physchem_filters", {})

    mask = df["p_active"] >= posterior_cutoff
    if require_hdi_positive:
        mask &= df["beta_hdi_lower"] > 0.0
    for col, (low, high) in physchem_filters.items():
        if col in df.columns:
            mask &= df[col].between(low, high)
    hits_df = df[mask].copy()

    # Output applied filters
    
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

    # Check if no hits, gracefully fail and warn
    if hits_df.empty:
        logger.warning("[del_hit_picker] No hits found after filtering.")
        
        # Return empty structures, but still in the same format
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


