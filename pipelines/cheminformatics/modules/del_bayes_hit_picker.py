import duckdb
import os

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
    description="Fit hierarchical Bayesian model to synthon counts per condition."
)
def del_bayesian_model(config: dict, data: dict) -> dict:
    params = config.get("del_bayesian_model", {})
    df = data.get("df")

    if df is None or len(df) == 0:
        parquet_path = data.get("input_file")
        if parquet_path and Path(parquet_path).exists():
            df = pd.read_parquet(parquet_path)
        else:
            raise ValueError("[del_bayesian_model] Input dataframe empty and no fallback found.")

    if "total_count" in df.columns:
        df = df.rename(columns={"total_count": "pre_count"})
        df["post_count"] = df["pre_count"]

    required_cols = {"NsynthonID", "condition", "pre_count", "post_count"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    # Encode categorical
    syn_codes, syn_idx = np.unique(df["NsynthonID"], return_inverse=True)
    cond_codes, cond_idx = np.unique(df["condition"], return_inverse=True)

    pre = df["pre_count"].values.astype(float)
    post = df["post_count"].values.astype(int)

    # Spike test
    spike_cfg = params.get("spike_test", {})
    if spike_cfg.get("enable", False):
        rng = np.random.default_rng(params.get("random_seed", 42))
        spike_idx = rng.choice(len(df), size=min(spike_cfg.get("synthons", 10), len(df)), replace=False)
        df.loc[spike_idx, "post_count"] = df.loc[spike_idx, "pre_count"] * spike_cfg.get("multiplier", 10)
        post = df["post_count"].values.astype(int)

    overdisp = params.get("overdispersion", 1.5)
    draws = params.get("draws", 1000)
    tune = params.get("tune", 1000)
    target_accept = params.get("target_accept", 0.9)
    random_seed = params.get("random_seed", 42)

    # Build model
    with pm.Model() as model:
        alpha_cond = pm.Normal("alpha_cond", mu=0.0, sigma=1.0, shape=len(cond_codes))
        sigma_syn = pm.HalfNormal("sigma_syn", 1.0)
        beta_syn = pm.Normal("beta_syn", mu=0.0, sigma=sigma_syn, shape=len(syn_codes))

        log_mu = alpha_cond[cond_idx] + beta_syn[syn_idx] + np.log(pre)
        mu = pm.math.exp(log_mu)

        pm.NegativeBinomial("obs", mu=mu, alpha=overdisp, observed=post)

        trace = pm.sample(
            draws=draws,
            tune=tune,
            target_accept=target_accept,
            random_seed=random_seed,
            chains=4,
            progressbar=True,
        )

    beta_post = trace.posterior["beta_syn"]
    p_active = (beta_post > 0).mean(dim=("chain", "draw")).values

    # HDI
    hdi = az.hdi(beta_post, hdi_prob=0.95)
    hdi_da = hdi["beta_syn"]
    hdi_vals = hdi_da.values
    hdi_axis = hdi_vals.shape.index(2)
    beta_hdi_lower = np.take(hdi_vals, 0, axis=hdi_axis)
    beta_hdi_upper = np.take(hdi_vals, 1, axis=hdi_axis)
    beta_mean = beta_post.mean(dim=("chain", "draw")).values

    # Build summary per (synthon, condition)
    summary_df = pd.DataFrame({
        "NsynthonID": df["NsynthonID"],
        "condition": df["condition"],
        "beta_mean": beta_mean[syn_idx],
        "beta_hdi_lower": beta_hdi_lower[syn_idx],
        "beta_hdi_upper": beta_hdi_upper[syn_idx],
        "p_active": p_active[syn_idx],
    })

    return {"df": summary_df, "trace": trace}

@register_task(
    "del_hit_picker",
    category="DEL",
    description="Select DEL synthons based on Bayesian posteriors."
)
def del_hit_picker(config: dict, data: dict) -> dict:

    import os

    params = config.get("del_hit_picker", {})
    df = data.get("df")

    if df is None or len(df) == 0:
        raise ValueError("[del_hit_picker] Input dataframe is empty.")

    # Convert to Polars
    if isinstance(df, pd.DataFrame):
        pl_df = pl.from_pandas(df)
    elif isinstance(df, pl.DataFrame):
        pl_df = df
    else:
        raise TypeError(f"Unsupported df type: {type(df)}")

    # Required columns
    required = {"NsynthonID", "beta_mean", "beta_hdi_lower", "beta_hdi_upper", "p_active", "condition"}
    missing = required - set(pl_df.columns)
    if missing:
        raise ValueError(f"[del_hit_picker] Missing columns: {missing}")

    posterior_cutoff = params.get("posterior_cutoff", 0.95)
    require_hdi_positive = params.get("require_hdi_positive", True)
    physchem_filters = params.get("physchem_filters", {})

    # Logging summaries
    logger.info("[del_hit_picker] Total synthons received: %d", pl_df.height)
    logger.info(
        "[del_hit_picker] p_active summary: min=%.3f median=%.3f max=%.3f",
        pl_df["p_active"].min(),
        pl_df["p_active"].median(),
        pl_df["p_active"].max()
    )
    logger.info(
        "[del_hit_picker] beta_hdi_lower summary: min=%.3f median=%.3f max=%.3f",
        pl_df["beta_hdi_lower"].min(),
        pl_df["beta_hdi_lower"].median(),
        pl_df["beta_hdi_lower"].max()
    )

    # Apply Bayesian filters
    mask = pl_df["p_active"] >= posterior_cutoff
    if require_hdi_positive:
        mask &= pl_df["beta_hdi_lower"] > 0.0

    filtered_hits = pl_df.filter(mask)
    logger.info("[del_hit_picker] Synthons after posterior/HDI filter: %d", filtered_hits.height)

    # Apply physchem filters
    for col, (low, high) in physchem_filters.items():
        if col in filtered_hits.columns:
            filtered_hits = filtered_hits.filter((pl.col(col) >= low) & (pl.col(col) <= high))
    logger.info("[del_hit_picker] Synthons after physchem filters: %d", filtered_hits.height)

    # Check for zero hits and log possible reasons
    if filtered_hits.height == 0:
        reasons = []
        if (pl_df["p_active"] >= posterior_cutoff).sum() == 0:
            reasons.append(f"posterior_cutoff ({posterior_cutoff}) too strict")
        if require_hdi_positive and (pl_df["beta_hdi_lower"] > 0).sum() == 0:
            reasons.append("HDI lower bound > 0 removed all synthons")
        for col, (low, high) in physchem_filters.items():
            if col in pl_df.columns:
                phys_mask = (pl_df[col] >= low) & (pl_df[col] <= high)
                if phys_mask.sum() == 0:
                    reasons.append(f"physchem filter {col} [{low}, {high}] removed all synthons")
        if not reasons:
            reasons.append("unknown, possibly small dataset or no real hits")
        logger.warning("[del_hit_picker] ZERO HITS Possible reasons: %s", "; ".join(reasons))

    # Sorting full hits
    sort_by = params.get("sort_by", ["p_active", "beta_mean"])
    ascending = params.get("ascending", [False, False])
    if len(sort_by) != len(ascending):
        raise ValueError("[del_hit_picker] 'sort_by' and 'ascending' must have same length.")
    descending = [not a for a in ascending]
    hits = filtered_hits.sort(by=sort_by, descending=descending)
    logger.info("[del_hit_picker] Hits after sorting: %d", hits.height)

    # Top-hit-per-synthon summary
    top_hits = (
        hits
        .sort(["p_active", "beta_mean"], descending=[True, True])
        .group_by("NsynthonID")
        .agg([
            pl.first("condition").alias("top_condition"),
            pl.first("beta_mean").alias("top_beta_mean"),
            pl.first("beta_hdi_lower").alias("top_beta_hdi_lower"),
            pl.first("beta_hdi_upper").alias("top_beta_hdi_upper"),
            pl.first("p_active").alias("top_p_active")
        ])
        .sort("top_p_active", descending=True)
    )
    logger.info("[del_hit_picker] Top-hit-per-synthon summary rows: %d", top_hits.height)

    # --- Save CSVs ---
    output_cfg = params.get("output", {})
    out_dir = Path(output_cfg.get("directory", "outputs/del_hits"))
    out_dir.mkdir(parents=True, exist_ok=True)

    full_csv_path = out_dir / output_cfg.get("filename", "wuxi_del_hits.csv")
    top_csv_path = out_dir / output_cfg.get("top_hits_filename", "wuxi_del_top_hits.csv")

    hits.to_pandas().to_csv(full_csv_path, index=False)
    top_hits.to_pandas().to_csv(top_csv_path, index=False)

    logger.info("[del_hit_picker] Full hits CSV saved: %s", full_csv_path)
    logger.info("[del_hit_picker] Top-hit summary CSV saved: %s", top_csv_path)

    return {
        "df": hits.to_pandas(),
        "top_hits_df": top_hits.to_pandas()
    }



