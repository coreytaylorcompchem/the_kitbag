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
from modules.utils.plot_del_hits import plot_del_hits

logger = setup_logger(__name__, debug_mode=False, simple_format=True)


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
    description="Select DEL synthons based on Bayesian posteriors and generate plots."
)
def del_hit_picker(config: dict, data: dict) -> dict:
    import warnings
    from pathlib import Path
    from modules.utils.plot_del_hits import plot_del_hits
    import polars as pl
    import pandas as pd

    params = config.get("del_hit_picker", {})
    df = data.get("df")  # Bayesian posterior output
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
            if c not in {"NsynthonID", "condition", "total_count", "n_compounds"} 
            and pd.api.types.is_numeric_dtype(physchem_df[c])
        ]

        if physchem_cols:
            merge_cols = ["NsynthonID", "condition"] + physchem_cols
            df = df.merge(physchem_df[merge_cols], on=["NsynthonID", "condition"], how="left")

            # Remove any physchem columns that became entirely NaN
            physchem_cols = [c for c in physchem_cols if c in df.columns and df[c].notna().any()]

    # ---------------- Apply Bayesian and physchem filters ----------------
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

    # ---------------- Sorting ----------------
    sort_by = params.get("sort_by", ["p_active", "beta_mean"])
    ascending = params.get("ascending", [False, False])
    hits_df = hits_df.sort_values(by=sort_by, ascending=ascending)

    # ---------------- Top-hit-per-synthon summary ----------------
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

    # ---------------- Save CSVs ----------------
    out_dir = Path(params.get("output", {}).get("directory", "outputs/del_hits"))
    out_dir.mkdir(parents=True, exist_ok=True)
    full_csv_path = out_dir / params.get("output", {}).get("filename", "wuxi_del_hits.csv")
    top_csv_path = out_dir / params.get("output", {}).get("top_hits_filename", "wuxi_del_top_hits.csv")

    hits_df.to_csv(full_csv_path, index=False)
    top_hits_df.to_csv(top_csv_path, index=False)

    logger.info("[del_hit_picker] Full hits CSV saved: %s", str(full_csv_path))
    logger.info("[del_hit_picker] Top-hit summary CSV saved: %s", str(top_csv_path))

    # ---------------- Generate plots ----------------
    plot_dir = params.get("output", {}).get("directory", "outputs/plots")
    Path(plot_dir).mkdir(parents=True, exist_ok=True)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=FutureWarning)

        # Only pass numeric physchem columns that exist in hits_df
        valid_physchem_cols = [c for c in physchem_cols if c in hits_df.columns and pd.api.types.is_numeric_dtype(hits_df[c])]

        plot_files = plot_del_hits(
            hits_df,
            output_dir=plot_dir,
            top_n=params.get("top_n_per_condition", 10),
            physchem_cols=valid_physchem_cols
        )

    # Convert PosixPath to str in logs
    plot_files_str = {k: str(v) for k, v in plot_files.items()}
    logger.info("[del_hit_picker] Plots generated: %s", plot_files_str)

    return {
        "df": hits_df,
        "top_hits_df": top_hits_df,
        "plot_files": plot_files_str
    }


