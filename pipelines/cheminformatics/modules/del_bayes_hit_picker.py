import duckdb

import polars as pl
import pandas as pd
import numpy as np
import pymc as pm
import arviz as az
from pathlib import Path

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)


@register_task(
    "del_stream_aggregate_counts",
    category="DEL",
    description="Stream DEL CSV chunks and aggregate synthon-level counts using DuckDB."
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

    agg_exprs = [f"SUM({count_col}) AS total_count", "COUNT(*) AS n_compounds"]
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

@register_task(
    "del_bayesian_model",
    category="DEL",
    description="Fit hierarchical Bayesian model to synthon counts."
)
def del_bayesian_model(config: dict, data: dict) -> dict:
    params = config.get("del_bayesian_model", {})
    
    df = data.get("df")

    # Fallback to Parquet if df is None or empty
    if df is None or len(df) == 0:
        parquet_path = data.get("input_file")  # <-- use input_file
        if parquet_path and Path(parquet_path).exists():
            df = pd.read_parquet(parquet_path)
        else:
            raise ValueError("[del_bayesian_model] Input dataframe is empty and no Parquet path available.")

    # Rename total_count -> pre_count, post_count same as pre_count
    if "total_count" in df.columns:
        df = df.rename(columns={"total_count": "pre_count"})
        df["post_count"] = df["pre_count"]

    required_cols = {"NsynthonID", "condition", "pre_count", "post_count"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Encode categorical variables
    synthon_codes, synthon_idx = np.unique(df["NsynthonID"], return_inverse=True)
    cond_codes, cond_idx = np.unique(df["condition"], return_inverse=True)

    pre = df["pre_count"].values.astype(float)
    post = df["post_count"].values.astype(int)

    if np.any(pre <= 0):
        raise ValueError("Found pre_count <= 0; filter these rows upstream.")

    # Model hyperparameters
    overdisp = params.get("overdispersion", 1.5)
    draws = params.get("draws", 1000)
    tune = params.get("tune", 1000)
    target_accept = params.get("target_accept", 0.9)
    random_seed = params.get("random_seed", 42)

    logger.info(f"[del_bayesian_model] Fitting model ({len(synthon_codes)} synthons, {len(cond_codes)} conditions)")

    with pm.Model() as model:
        alpha_cond = pm.Normal("alpha_cond", mu=0.0, sigma=1.0, shape=len(cond_codes))
        sigma_syn = pm.HalfNormal("sigma_syn", 1.0)
        beta_syn = pm.Normal("beta_syn", mu=0.0, sigma=sigma_syn, shape=len(synthon_codes))
        log_mu = alpha_cond[cond_idx] + beta_syn[synthon_idx] + np.log(pre)
        mu = pm.math.exp(log_mu)
        pm.NegativeBinomial("obs", mu=mu, alpha=overdisp, observed=post)

        trace = pm.sample(
            draws=draws,
            tune=tune,
            target_accept=target_accept,
            random_seed=random_seed,
            chains=4,
            progressbar=True
        )

    beta_post = trace.posterior["beta_syn"]
    p_active = (beta_post > 0).mean(dim=("chain", "draw")).values
    hdi = az.hdi(beta_post, hdi_prob=0.95).to_array().values

    summary_df = pd.DataFrame({
        "NsynthonID": synthon_codes,
        "beta_mean": beta_post.mean(dim=("chain", "draw")).values,
        "beta_hdi_lower": hdi[0],
        "beta_hdi_upper": hdi[1],
        "p_active": p_active
    })

    return {"df": summary_df, "trace": trace}


@register_task(
    "del_hit_picker",
    category="DEL",
    description="Select DEL synthons based on Bayesian posteriors."
)
def del_hit_picker(config: dict, data: dict) -> dict:

    params = config.get("del_hit_picker", {})
    df = data.get("df")

    if df is None or len(df) == 0:
        raise ValueError("[del_hit_picker] Input dataframe is empty.")

    # Convert to Polars for safe filtering
    if isinstance(df, pd.DataFrame):
        pl_df = pl.from_pandas(df)
    elif isinstance(df, pl.DataFrame):
        pl_df = df
    else:
        raise TypeError(f"Unsupported df type: {type(df)}")

    # Required columns
    required = {"NsynthonID", "beta_mean", "beta_hdi_lower", "beta_hdi_upper", "p_active"}
    missing = required - set(pl_df.columns)
    if missing:
        raise ValueError(f"[del_hit_picker] Missing columns: {missing}")

    posterior_cutoff = params.get("posterior_cutoff", 0.95)
    require_hdi_positive = params.get("require_hdi_positive", True)

    # Mask
    mask = pl_df["p_active"] >= posterior_cutoff
    if require_hdi_positive:
        mask &= pl_df["beta_hdi_lower"] > 0.0

    hits = pl_df.filter(mask)

    # Optional physchem filtering
    physchem_filters = params.get("physchem_filters", {})
    for col, (low, high) in physchem_filters.items():
        if col in hits.columns:
            hits = hits.filter((pl.col(col) >= low) & (pl.col(col) <= high))

    # Sorting
    sort_by = params.get("sort_by", ["p_active", "beta_mean"])
    ascending = params.get("ascending", [False, False])
    hits = hits.sort(by=sort_by, reverse=[not a for a in ascending])

    return {"df": hits.to_pandas()}