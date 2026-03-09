import os
import psutil
import time
from rdkit import Chem

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

# ----------------------------
# IMPORTS AFTER THREAD LIMITS
# ----------------------------
import duckdb
import warnings

from collections import Counter
from math import log
import polars as pl
import pandas as pd
import numpy as np
import pymc as pm
import arviz as az
from pathlib import Path
from itertools import combinations

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger
from modules.utils.plot_del_hits import plot_del_hits, reduce_hits_for_plotting
from modules.utils.plot_del_hits import plot_del_qc_metrics
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
    
    # Propagate product (SMILES) through Bayes model
    if "product" in df.columns:
        agg_exprs.append("ANY_VALUE(product) AS product")
        
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
    description="Fit hierarchical Bayesian model to synthon enrichment per condition using ADVI with low-count filtering."
)
def del_bayesian_model(config: dict, data: dict) -> dict:

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
        
        # HIERARCHICAL PARTIAL POOLING 
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

    # ----------------------------
    # Propagate product SMILES
    # ----------------------------
    if "product" in df.columns:
        smiles_df = (
            df.groupby(["NsynthonID", "condition"], as_index=False)
            .agg(product=("product", "first"))
        )

        summary_df = summary_df.merge(
            smiles_df,
            on=["NsynthonID", "condition"],
            how="left"
        )

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


    # Propagate enrichment for QC plotting

    summary_df["enrichment_mean"] = (
        model_df.set_index(["NsynthonID", "condition"])["enrichment_mean"]
        .reindex(summary_df.set_index(["NsynthonID", "condition"]).index)
        .values
    )

    return {"df": summary_df, "trace": trace, "physchem_cols": physchem_cols}

@register_task(
    "del_hit_picker",
    category="DEL",
    description="Select DEL synthons based on Bayesian posteriors and generate plots."
)
def del_hit_picker(config: dict, data: dict) -> dict:
    
    params = config.get("del_hit_picker", {})
    df = data.get("df")  # Output from del_bayesian_model

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
        physchem_cols = [c for c in physchem_cols if c in df.columns]
    
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

    logger.info(
        "[del_hit_picker] Physchem filters requested (YAML keys): %s",
        list(physchem_filters.keys())
    )

    if physchem_filters:
        logger.info(
            "[del_hit_picker] Physchem filters requested from YAML: %s",
            physchem_filters
        )
    else:
        logger.info("[del_hit_picker] No physchem filters provided in YAML")

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

    for base_col, (low, high) in physchem_filters.items():

        # Resolve column name
        if base_col in df.columns:
            col = base_col
        elif f"{base_col}_mean" in df.columns:
            col = f"{base_col}_mean"
        elif f"{base_col}_mean_mean" in df.columns:
            col = f"{base_col}_mean_mean"
        else:
            logger.warning(
                "[del_hit_picker] Physchem filter '%s' not found (checked %s, %s_mean, %s_mean_mean) — skipping",
                base_col, base_col, base_col, base_col
            )
            continue

        before = mask.sum()
        mask_phys = df[col].between(low, high)
        mask &= mask_phys
        after = mask.sum()

        logger.info(
            "[del_hit_picker] Physchem filter applied: %s → %s in [%.2f, %.2f] | %d → %d rows",
            base_col, col, low, high, before, after
        )

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

    # umap_features = [c for c in physchem_cols if c in hits_df.columns]

    # logger.info(
    #     "[del_hit_picker] UMAP features: %s",
    #     umap_features
    # )

    # hits_df = add_umap_clusters(
    #     hits_df,
    #     feature_cols=umap_features,
    #     n_neighbors=30,
    #     min_dist=0.15,
    #     cluster_method="hdbscan",
    # )

    # ----------------------------
    # UMAP embedding on unique synthons
    # ----------------------------

    # Aggregate to one row per NsynthonID
    umap_base_df = hits_df.groupby("NsynthonID", as_index=False).agg({
        "beta_mean": "max",
        "p_active": "max",
        **{c: "mean" for c in physchem_cols if c in hits_df.columns}
    })

    # Extract library
    umap_base_df["library"] = umap_base_df["NsynthonID"].str.split("-", n=1).str[0]

    # Compute library weights (inverse frequency)
    lib_counts = umap_base_df["library"].value_counts()
    umap_base_df["umap_weight"] = 1.0 / lib_counts[umap_base_df["library"]].values

    # Sample for UMAP (weighted, unique synthons only)
    n_umap_samples = min(len(umap_base_df), params.get("max_umap_points", 10000))
    umap_sample_df = umap_base_df.sample(
        n=n_umap_samples,
        weights="umap_weight",
        random_state=42
    ).reset_index(drop=True)

    logger.info(
        "[del_hit_picker] Using %d unique synthons for UMAP (%d libraries, proportional sampling)",
        len(umap_sample_df),
        umap_sample_df["library"].nunique()
    )

    # Compute UMAP embedding
    umap_features = ["beta_mean", "p_active"] + [c for c in physchem_cols if c in umap_sample_df.columns]

    umap_sample_df = add_umap_clusters(
        umap_sample_df,
        feature_cols=umap_features,
        n_neighbors=30,
        min_dist=0.15,
        cluster_method="hdbscan",
    )

    # Merge UMAP coordinates back safely into full hits_df
    # Only one row per NsynthonID to avoid dups
    hits_df = hits_df.merge(
        umap_sample_df[["NsynthonID", "umap_x", "umap_y", "cluster_id"]],
        on="NsynthonID",
        how="left"
    )

    # ----------------------------
    # Reduce hits for plotting only
    # ----------------------------
    
    max_plot_hits = params.get("max_plot_hits", 50000) # TODO: ADD THIS TO YAML

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
            hits_df=reduced_hits_df,
            full_hits_df=hits_df,
            output_dir=plot_dir,
            top_n=params.get("top_n_per_condition", 10),
            physchem_cols=physchem_cols,
            physchem_filters=physchem_filters,  # NEW
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
    "del_qc_metrics",
    category="DEL",
    description="Ad-hoc QC diagnostics for DEL Bayesian model."
)
def del_qc_metrics(config: dict, data: dict) -> dict:

    logger.info("[del_qc_metrics] Running ad-hoc DEL QC metrics")

    df = data.get("df")
    trace = data.get("trace")

    if df is None or trace is None:
        logger.warning("[del_qc_metrics] Missing df or trace - skipping QC")
        return {"qc_metrics": {}}

    df = df.copy()

    qc = {}

    # ----------------------------
    # 1. Global signal sanity
    # ----------------------------
    beta = df["beta_mean"].values
    p_active = df["p_active"].values

    qc["global"] = {
        "beta_mean_mean": float(np.mean(beta)),
        "beta_mean_std": float(np.std(beta)),
        "frac_p_active_gt_0.5": float((p_active > 0.5).mean()),
        "frac_p_active_gt_0.95": float((p_active > 0.95).mean()),
        "n_rows": int(len(df)),
    }

    if qc["global"]["frac_p_active_gt_0.95"] < 0.01:
        logger.warning(
            "[QC] Very few strong hits (p_active > 0.95: %.2f%%)",
            100 * qc["global"]["frac_p_active_gt_0.95"]
        )

    # ----------------------------
    # 2. Posterior uncertainty
    # ----------------------------
    hdi_width = df["beta_hdi_upper"] - df["beta_hdi_lower"]

    qc["posterior_uncertainty"] = {
        "median_hdi_width": float(np.median(hdi_width)),
        "p90_hdi_width": float(np.quantile(hdi_width, 0.9)),
        "frac_p_active_but_hdi_crosses_0": float(
            ((df["p_active"] > 0.95) & (df["beta_hdi_lower"] <= 0)).mean()
        ),
    }

    if qc["posterior_uncertainty"]["frac_p_active_but_hdi_crosses_0"] > 0.3:
        logger.warning(
            "[QC] Many high p_active hits have HDIs crossing zero (%.1f%%)",
            100 * qc["posterior_uncertainty"]["frac_p_active_but_hdi_crosses_0"]
        )

    # ----------------------------
    # 3. Condition-level diagnostics
    # ----------------------------
    cond_stats = (
        df.groupby("condition")
        .agg(
            mean_beta=("beta_mean", "mean"),
            std_beta=("beta_mean", "std"),
            frac_active=("p_active", lambda x: (x > 0.95).mean()),
            n=("beta_mean", "size"),
        )
        .reset_index()
    )

    qc["conditions"] = cond_stats.to_dict(orient="records")

    if cond_stats["frac_active"].max() > 3 * cond_stats["frac_active"].median():
        logger.warning(
            "[QC] One condition dominates hit rate — possible condition-specific artifact"
        )

    # ----------------------------
    # 4. Library dominance
    # ----------------------------
    if "NsynthonID" in df.columns:
        libs = df["NsynthonID"].str.split("-", n=1).str[0]
        hit_libs = libs[p_active > 0.95]

        lib_counts = Counter(hit_libs)
        total = sum(lib_counts.values())

        def shannon_entropy(counter):
            return -sum((c/total) * log(c/total) for c in counter.values() if c > 0)

        qc["libraries"] = {
            "n_libraries_with_hits": len(lib_counts),
            "top_library_frac": (
                max(lib_counts.values()) / total if total > 0 else 0.0
            ),
            "library_entropy": shannon_entropy(lib_counts) if total > 0 else 0.0,
        }

        if qc["libraries"]["top_library_frac"] > 0.6:
            logger.warning(
                "[QC] One library contributes %.1f%% of hits — possible bias",
                100 * qc["libraries"]["top_library_frac"]
            )

    # ----------------------------
    # 5. Variance decomposition
    # ----------------------------
    try:
        var_summary = az.summary(
            trace,
            var_names=["sigma_obs", "sigma_cond"],
            kind="stats"
        )

        sigma_obs = var_summary.loc["sigma_obs", "mean"]
        sigma_cond = var_summary.loc["sigma_cond", "mean"]

        qc["variance"] = {
            "sigma_obs_mean": float(sigma_obs),
            "sigma_cond_mean": float(sigma_cond),
            "sigma_cond_over_obs": float(sigma_cond / sigma_obs),
        }

        if sigma_obs > 2 * sigma_cond:
            logger.warning(
                "[QC] Observation noise dominates condition effects (sigma_obs >> sigma_cond)"
            )

    except Exception as e:
        logger.warning("[QC] Could not compute variance diagnostics: %s", e)

    # ----------------------------
    # 6. QC plots
    # ----------------------------
    try:
        plot_dir = (
            config
            .get("output", {})
            .get("plots_directory", "outputs/plots_qc")
        )

        plot_del_qc_metrics(
            df,
            output_dir=plot_dir
        )

        logger.info("[del_qc_metrics] QC plots written to %s", plot_dir)

    except Exception as e:
        logger.warning(
            "[del_qc_metrics] QC plotting failed (non-fatal): %s", e
        )


    logger.info("[del_qc_metrics] QC metrics completed")

    return {"qc_metrics": qc, "df": df, "trace": trace}

@register_task(
    "del_smarts_filter",
    category="DEL",
    description="Match molecules against toxic/reactive SMARTS patterns, assign traffic-light severity."
)
def del_smarts_filter(config: dict, data: dict) -> dict:

    params = config.get("del_smarts_filter", {})

    smarts_file = params.get("smarts_file")
    smiles_col = params.get("smiles_column", "product")

    severity_col = params.get("output_columns", {}).get("severity", "smarts_severity")
    matches_col = params.get("output_columns", {}).get("matches", "smarts_matches")

    df = data.get("df")

    if df is None or len(df) == 0:
        logger.warning("[del_smarts_filter] Empty dataframe")
        return {"df": df}

    df = df.copy()

    if smiles_col not in df.columns:
        logger.warning("[del_smarts_filter] SMILES column '%s' not found — skipping SMARTS filter", smiles_col)
        df[severity_col] = "GREEN"
        df[matches_col] = ""
        return {"df": df}

    logger.info("[del_smarts_filter] Loading SMARTS from %s", smarts_file)

    smarts_df = pd.read_csv(smarts_file)

    patterns = []

    for _, row in smarts_df.iterrows():
        patt = Chem.MolFromSmarts(row["SMARTS"])
        if patt is None:
            logger.warning("[del_smarts_filter] Invalid SMARTS skipped: %s", row["SMARTS"])
            continue

        patterns.append({
            "name": row["PATTERN_NAME"],
            "pattern": patt,
            "severity": row["SEVERITY_COMMENT"].upper()
        })

    logger.info("[del_smarts_filter] Loaded %d SMARTS patterns", len(patterns))

    severities = []
    matches_list = []

    for smi in df[smiles_col].astype(str):

        mol = Chem.MolFromSmiles(smi)

        if mol is None:
            severities.append("GREEN")
            matches_list.append("")
            continue

        matched_names = []
        matched_severities = []

        for p in patterns:
            if mol.HasSubstructMatch(p["pattern"]):
                matched_names.append(p["name"])
                matched_severities.append(p["severity"])

        if "REJECT" in matched_severities:
            sev = "RED"
        elif "FLAG" in matched_severities:
            sev = "ORANGE"
        else:
            sev = "GREEN"

        severities.append(sev)
        matches_list.append(",".join(matched_names))

    df[severity_col] = severities
    df[matches_col] = matches_list

    logger.info(
        "[del_smarts_filter] Results: GREEN=%d ORANGE=%d RED=%d",
        (df[severity_col] == "GREEN").sum(),
        (df[severity_col] == "ORANGE").sum(),
        (df[severity_col] == "RED").sum(),
    )

    return {"df": df}

@register_task(
    "del_compare_datasets",
    category="DEL",
    description="Compare raw DEL data to processed DEL data and report missing synthons."
)
def del_compare_datasets(config: dict, data: dict) -> dict:
    params = config.get("del_compare_datasets", {})
    processed_input = params.get("processed_input")
    check_values = params.get("check_values", False)
    output_cfg = params.get("output", {})

    out_dir = Path(output_cfg.get("directory", "outputs/del_compare"))
    out_dir.mkdir(parents=True, exist_ok=True)

    # -------------------
    # Load RAW (long format)
    # -------------------
    input_file = data.get("input_file")
    if not input_file:
        logger.error("No input_file found in postprocess data.")
        return {}

    raw_df = pd.read_parquet(input_file)
    if raw_df.empty:
        logger.warning("Aggregated raw dataframe is empty.")
        return {}

    raw_df["NsynthonID"] = raw_df["NsynthonID"].astype(str).str.strip()

    # Ensure counts are numeric
    if "post_count" in raw_df.columns:
        # Ensure counts numeric
        raw_df["post_count"] = pd.to_numeric(raw_df["post_count"], errors="coerce").fillna(0)

        # Pivot counts wide
        raw_count_wide = raw_df.pivot_table(
            index="NsynthonID",
            columns="condition",
            values="post_count",
            aggfunc="sum",
            fill_value=0
        )

        raw_count_wide.columns.name = None

        # -------------------
        # Filter 1: remove compounds with zero total counts
        # -------------------
        total_counts = raw_count_wide.sum(axis=1)
        mask_nonzero_total = total_counts > 0

        logger.info(f"Filtered synthons with count = 0 in all conditions."
            f"Remaining after threshold filtering: {len(mask_nonzero_total)}"
            )

        # -------------------
        # Filter 2: remove compounds with < 2 in ALL conditions
        # -------------------
        threshold = 2
        mask_threshold = (raw_count_wide >= threshold).any(axis=1)

        # Combine filters
        valid_ids = raw_count_wide.index[mask_nonzero_total & mask_threshold]

        raw_df = raw_df[raw_df["NsynthonID"].isin(valid_ids)]

        logger.info(f"Filtered synthons with count < {threshold} in all conditions."
            f"Remaining after threshold filtering: {len(valid_ids)}"
            )

    raw_unique = raw_df[["NsynthonID"]].drop_duplicates()
    logger.info(f"Raw unique synthons: {len(raw_unique)}")

    # -------------------
    # Load processed (wide format)
    # -------------------
    processed_files = list(Path(processed_input).glob("*.csv"))
    if not processed_files:
        raise ValueError("No processed CSV files found.")

    processed_df = pd.concat(
        [pd.read_csv(f, low_memory=False) for f in processed_files],
        ignore_index=True
    )
    processed_df["NsynthonID"] = processed_df["NsynthonID"].astype(str).str.strip()
    processed_unique = processed_df[["NsynthonID"]].drop_duplicates()
    logger.info(f"Processed unique synthons: {len(processed_unique)}")

    # -------------------
    # Check for duplicates
    # -------------------
    raw_dupes = raw_df["NsynthonID"].duplicated().sum()
    processed_dupes = processed_df["NsynthonID"].duplicated().sum()
    if raw_dupes:
        logger.warning(f"Raw dataset has {raw_dupes} duplicate NsynthonIDs")
    if processed_dupes:
        logger.warning(f"Processed dataset has {processed_dupes} duplicate NsynthonIDs")

    # -------------------
    # Compare IDs
    # -------------------
    raw_ids = set(raw_unique["NsynthonID"])
    processed_ids = set(processed_unique["NsynthonID"])
    intersection = raw_ids & processed_ids
    missing_ids = raw_ids - processed_ids

    logger.info(f"Raw ID count: {len(raw_ids)}")
    logger.info(f"Processed ID count: {len(processed_ids)}")
    logger.info(f"Intersection size: {len(intersection)}")
    logger.info(f"Missing synthons: {len(missing_ids)}")

    # -------------------
    # Output missing IDs with per-condition data
    # -------------------
    if missing_ids:

        # Ensure enrichment numeric
        raw_df["enrichment_mean"] = pd.to_numeric(raw_df["enrichment_mean"], errors="coerce").fillna(0)

        # Pivot counts
        raw_count_wide = raw_df.pivot_table(
            index="NsynthonID",
            columns="condition",
            values="post_count",
            aggfunc="sum",
            fill_value=0
        )

        raw_count_wide.columns = [f"count_{c}" for c in raw_count_wide.columns]
        raw_count_wide.reset_index(inplace=True)

        # Pivot enrichment
        raw_enrich_wide = raw_df.pivot_table(
            index="NsynthonID",
            columns="condition",
            values="enrichment_mean",
            aggfunc="mean",
            fill_value=0
        )

        raw_enrich_wide.columns = [f"enrichment_{c}" for c in raw_enrich_wide.columns]
        raw_enrich_wide.reset_index(inplace=True)

        # Merge count + enrichment
        raw_wide = raw_count_wide.merge(raw_enrich_wide, on="NsynthonID", how="left")

        # merge product SMILES
        if "product" in raw_df.columns:
            product_df = (
                raw_df[["NsynthonID", "product"]]
                .drop_duplicates()
            )
            raw_wide = raw_wide.merge(product_df, on="NsynthonID", how="left")

        # Filter to missing
        missing_df = raw_wide[raw_wide["NsynthonID"].isin(missing_ids)].copy()

        missing_outfile = out_dir / output_cfg.get("filename_missing", "missing_in_processed.csv")
        missing_df.to_csv(missing_outfile, index=False)

        logger.info("Missing synthons exported with per-condition counts and enrichment.")

    else:
        logger.info("No missing synthons to export.")

    # -------------------
    # Value comparison (wide, condition-by-condition) - WIP
    # -------------------
    if check_values:
        raw_df.columns = raw_df.columns.str.strip()
        processed_df.columns = processed_df.columns.str.strip()

        # Identify copy/enrichment columns in processed wide dataset
        copy_cols = [c for c in processed_df.columns if c.startswith("copy(")]
        enrich_cols = [c for c in processed_df.columns if c.startswith("enrichment(")]

        # Map raw columns correctly
        raw_copy_col = "post_count"
        raw_enrich_col = "enrichment_mean"

        # Ensure raw columns are numeric
        if raw_copy_col in raw_df.columns:
            raw_df[raw_copy_col] = pd.to_numeric(raw_df[raw_copy_col], errors="coerce").fillna(0)
        if raw_enrich_col in raw_df.columns:
            raw_df[raw_enrich_col] = pd.to_numeric(raw_df[raw_enrich_col], errors="coerce").fillna(0)

        # Ensure processed columns are numeric
        for col in copy_cols + enrich_cols:
            processed_df[col] = pd.to_numeric(processed_df[col], errors="coerce").fillna(0)

        # Pivot raw data long -> wide
        if raw_copy_col in raw_df.columns:
            raw_copy_wide = raw_df.pivot_table(
                index="NsynthonID",
                columns="condition",
                values=raw_copy_col,
                aggfunc="sum",
                fill_value=0
            )
            raw_copy_wide.columns.name = None
            raw_copy_wide = raw_copy_wide.reset_index()
        else:
            logger.warning(f"Raw dataframe has no '{raw_copy_col}' column. Skipping copy comparison.")
            raw_copy_wide = pd.DataFrame({"NsynthonID": raw_df["NsynthonID"].unique()})

        if raw_enrich_col in raw_df.columns:
            raw_enrich_wide = raw_df.pivot_table(
                index="NsynthonID",
                columns="condition",
                values=raw_enrich_col,
                aggfunc="mean",
                fill_value=0
            )
            raw_enrich_wide.columns.name = None
            raw_enrich_wide = raw_enrich_wide.reset_index()
        else:
            logger.warning(f"Raw dataframe has no '{raw_enrich_col}' column. Skipping enrichment comparison.")
            raw_enrich_wide = pd.DataFrame({"NsynthonID": raw_df["NsynthonID"].unique()})

        # Map processed columns to raw condition names
        copy_col_map = {col.split("(")[-1][:-1].split("_", 1)[-1]: col for col in copy_cols}
        enrich_col_map = {col.split("(")[-1][:-1].split("_", 1)[-1]: col for col in enrich_cols}

        # Merge processed and raw pivoted data
        merged = processed_df[["NsynthonID"] + copy_cols + enrich_cols].merge(
            raw_copy_wide, on="NsynthonID", how="inner"
        )

        diff_records = []
        # Compute differences for copy
        for cond, proc_col in copy_col_map.items():
            if cond in raw_copy_wide.columns:
                merged[f"{proc_col}_diff"] = merged[proc_col].astype(float) - merged[cond].astype(float)
                diff_records.append(f"{proc_col}_diff")
        # Compute differences for enrichment
        for cond, proc_col in enrich_col_map.items():
            if cond in raw_enrich_wide.columns:
                merged[f"{proc_col}_diff"] = merged[proc_col].astype(float) - merged[cond].astype(float)
                diff_records.append(f"{proc_col}_diff")

        mismatch_outfile = out_dir / output_cfg.get("filename_mismatch", "value_mismatches.csv")
        merged[["NsynthonID"] + diff_records].to_csv(mismatch_outfile, index=False)
        logger.info(f"[COMPARE] Value mismatches exported per condition")

    return {}