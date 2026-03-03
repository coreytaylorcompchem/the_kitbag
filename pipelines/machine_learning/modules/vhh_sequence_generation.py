import random
from collections import defaultdict
from pathlib import Path
import joblib

import numpy as np
import pandas as pd
import torch
import time

from tqdm import tqdm

from itertools import combinations
import matplotlib.cm as cm

import esm
from catboost import CatBoostRegressor

from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

import umap
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm

from pipeline.task_registry import register_task
from modules.utils.plotting import (
    compute_multi_condition_pareto,
    plot_round_radar,
    generate_summary_plots,
    plot_learning_curves_per_property,
)
from modules.utils.settings_logs import log_pipeline_config
from modules.utils.eval_regression import save_evaluation, evaluate_with_bootstrap
from modules.utils.eval_roundwise import run_roundwise_evaluation

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def debug_numeric_array(name, arr, max_examples=5):
    logger.debug(f"[DEBUG] {name}: dtype={arr.dtype}, shape={arr.shape}")

    # Check if object dtype
    if arr.dtype == object:
        logger.error(f"[DEBUG] {name} is OBJECT dtype")

        # Find first few non-numeric entries
        bad = []
        for i, v in enumerate(arr):
            if not isinstance(v, (int, float, np.integer, np.floating)) and not pd.isna(v):
                bad.append((i, type(v), v))
            if len(bad) >= max_examples:
                break

        logger.error(f"[DEBUG] {name} sample non-numeric entries: {bad}")

    # Test np.isfinite safely
    try:
        mask = np.isfinite(arr)
        logger.debug(f"[DEBUG] {name} np.isfinite OK, finite_count={mask.sum()}")
    except Exception as e:
        logger.error(f"[DEBUG] np.isfinite FAILED on {name}: {e}")

@register_task("load_seed_dataset", category="VHH generation", description="Load seed sequences and initialize context")
def load_seed_dataset(config, context):
    """
    Load seed sequences and initialize context.

    Args:
        config (dict): Task-specific YAML slice.
        context (dict): Runtime context.
        full_config (dict, optional): Entire YAML config for logging workflow-level info.
    """

    # Pass the full YAML to the logger

    full_config = context.get("full_config", {})
    log_pipeline_config(config, context, full_config=full_config)

    # read in exp data csv

    df = pd.read_csv(config["csv_path"])
    
    # Force numeric only on the physchem columns
    df[config["multi_condition_props"]] = df[config["multi_condition_props"]].apply(pd.to_numeric, errors="coerce").astype(float)
    expected_len = config.get("expected_len", 125)
    df = df[df["sequence"].str.len() == expected_len].reset_index(drop=True)

    # Convert multi-condition property columns to numeric
    for prop in config["multi_condition_props"]:
        df[prop] = pd.to_numeric(df[prop], errors="coerce")

    # Add ID number

    df["vhh_num"] = df["vhh_num"].astype(str)

    seeds = df["sequence"].tolist()
    seed_ids = list(range(len(df)))
    seq_len = len(seeds[0])

    # CDR positions from YAML
    cdr_regions = config.get("cdr_regions")
    if cdr_regions is None:
        raise ValueError("cdr_regions must be explicitly defined in YAML for VHHs")
    if len(cdr_regions) != 3:
        raise ValueError("Expected exactly 3 CDR regions for VHHs (CDR1–3)")

    # Framework positions from YAML
    framework_ranges = config.get("framework_positions")
    if not framework_ranges:
        logger.warning("No framework_positions found in YAML! All positions will be mutable.")
        framework_ranges = []

    # Keep the original ranges for plotting
    framework_ranges_tuples = [(start, end) for start, end in framework_ranges]

    # Flatten for mutable positions
    framework_positions_flat = []
    for start, end in framework_ranges_tuples:
        framework_positions_flat.extend(range(start, end))

    mutable_positions = []
    for start, end in cdr_regions.values():
        mutable_positions.extend(range(start, end))
    mutable_positions = sorted(set(mutable_positions))

    context.update({
        "df_seeds": df,
        "seeds": seeds,
        "seed_ids": seed_ids,
        "multi_condition_props": config["multi_condition_props"],
        "mutable_positions": mutable_positions,
        "cdr_regions": cdr_regions,
        "framework_positions": framework_ranges_tuples,
        "current_seeds": seeds.copy(),
        "current_seed_ids": seed_ids.copy(),
        "current_ancestor_ids": df["ancestor_id"].tolist(),
        "round_stats": [],
        "centroid_history": defaultdict(list)
    })

    logger.info(f"Loaded {len(seeds)} seed sequences")
    logger.info(f"Expected sequence length: {seq_len}")
    logger.info(f"Framework positions count: {len(framework_positions_flat)}")
    logger.info(f"Mutable positions count: {len(mutable_positions)}")

    return context


@register_task("load_esm_model", category="VHH generation", description="Load pretrained ESM model")
def load_esm_model(config, context):
    device = torch.device(config.get("device", "cuda" if torch.cuda.is_available() else "cpu"))
    model_name = config.get("pretrained_model", "esm1b_t33_650M_UR50S")
    esm_model, alphabet = getattr(esm.pretrained, model_name)()
    batch_converter = alphabet.get_batch_converter()
    esm_model.eval()
    esm_model = esm_model.to(device)

    context.update({
        "esm_model": esm_model,
        "alphabet": alphabet,
        "batch_converter": batch_converter,
        "device": device,
        "esm_cache": {},
        "batch_size": config.get("batch_size", 32)
    })
    logger.info(f"Loaded ESM model '{model_name}' on device {device}")
    return context

@register_task("active_learning_rounds", category="VHH generation", description="Run AL rounds with UMAP and radar plots")
def active_learning_rounds(config, context):

    n_rounds = config.get("n_rounds", 3)
    batch_size = config.get("batch_size", 32)
    samples_per_seed = config.get("samples_per_seed", 100)
    max_mutants = config.get("max_mutants", 3)
    top_k = config.get("top_candidates_per_round", 1000)
    min_per_ancestor = config.get("min_per_ancestor", 5)
    noise_scale = config.get("noise_scale", 1.0)

    multi_condition_props = context["multi_condition_props"]
    catboost_params = config.get("catboost_params", {})

    current_seeds = context["current_seeds"]
    current_seed_ids = context["current_seed_ids"]
    current_ancestor_ids = context["current_ancestor_ids"]
    df_labeled = context["df_seeds"][["sequence", "ancestor_id", "vhh_num"] + multi_condition_props].copy()
    mutable_positions = context["mutable_positions"]
    esm_model = context["esm_model"]
    batch_converter = context["batch_converter"]
    device = context["device"]
    esm_cache = context["esm_cache"]

    data_dir = Path(config.get("data_dir", "data_vhh"))
    plots_dir = Path(config.get("plots_dir", "plots_vhh"))
    data_dir.mkdir(exist_ok=True)
    plots_dir.mkdir(exist_ok=True)

    seq_len = len(current_seeds[0])
    aa_list = list("ACDEFGHIKLMNPQRSTVWY")
    aa_to_int = {aa: i for i, aa in enumerate(aa_list)}

    centroid_history = context["centroid_history"]

    # ------------------------
    # Helpers
    # ------------------------
    def generate_mutants(seq):
        mutants = set()
        for k in range(1, max_mutants+1):
            for positions in combinations(mutable_positions, k):
                new_seq = list(seq)
                for p in positions:
                    new_seq[p] = random.choice(aa_list)
                mutants.add(''.join(new_seq))
                if len(mutants) >= samples_per_seed:
                    break
            if len(mutants) >= samples_per_seed:
                break
        return list(mutants)

    @torch.no_grad()
    def esm_embed_cached(sequences):
        all_embeddings = []
        to_compute = []
        compute_idx = []

        for i, seq in enumerate(sequences):
            if seq in esm_cache:
                all_embeddings.append(esm_cache[seq])
            else:
                all_embeddings.append(None)
                to_compute.append(seq)
                compute_idx.append(i)

        for i in tqdm(range(0, len(to_compute), batch_size), desc=f"Embedding batch {round_idx}", leave=False):
            batch_seqs = to_compute[i:i+batch_size]
            data = [("seq", seq) for seq in batch_seqs]
            _, _, batch_tokens = batch_converter(data)
            batch_tokens = batch_tokens.to(device)
            with torch.cuda.amp.autocast():
                results = esm_model(batch_tokens, repr_layers=[33], return_contacts=False)
            token_reps = results["representations"][33]
            for j, seq in enumerate(batch_seqs):
                emb = token_reps[j, 1:len(seq)+1].mean(0).cpu().numpy()
                idx = compute_idx[i+j]
                all_embeddings[idx] = emb
                esm_cache[seq] = emb
            del batch_tokens, results, token_reps
            torch.cuda.empty_cache()

        return np.vstack(all_embeddings)

    # ------------------------
    # PCA setup
    # ------------------------
    pca = None

    for round_idx in range(1, n_rounds+1):
        round_start_time = time.perf_counter()
        logger.info(f"=== Round {round_idx} ===")

        # Embed labeled sequences
        X_train_full = esm_embed_cached(df_labeled["sequence"].tolist())

        # Ensure multi-condition props are numeric (only those columns)
        df_labeled[multi_condition_props] = df_labeled[multi_condition_props].apply(pd.to_numeric, errors="coerce")
        y_train = df_labeled[multi_condition_props].values

        logger.debug("[DEBUG] Checking whether y_train matrix is numeric")
        debug_numeric_array("y_train_flat", y_train.flatten())

        logger.info(f"Embedding {len(df_labeled)} labeled sequences...")

        if pca is None:
            n_samples, _ = X_train_full.shape
            pca_components = min(config.get("n_components", 512), n_samples)
            pca = PCA(n_components=pca_components, random_state=config.get("pca_seed", 42))
            pca.fit(X_train_full)
        X_train_pca = pca.transform(X_train_full)

        # ------------------------
        # Train CatBoost per property
        # ------------------------
        logger.info("Training CatBoost models for properties...")
        models = {}
        n_ensemble = config.get("n_model_ensemble", 5)      

        logger.info(f"Number of model ensembles: {n_ensemble}")

        for i, prop in enumerate(multi_condition_props):

            logger.info(f"Training model for: {prop}")

            y_i = y_train[:, i]

            # Debug numeric columns

            logger.debug(f"[DEBUG] Property {prop}")
            debug_numeric_array(f"{prop}_y_i_before_mask", y_i)
            
            # ensure float dtype
            if not np.issubdtype(y_i.dtype, np.floating):
                y_i = y_i.astype(float)
            try:
                mask = np.isfinite(y_i)
            except Exception as e:
                logger.error(f"[FATAL DEBUG] np.isfinite crashed for property {prop}")
                debug_numeric_array(f"{prop}_y_i_CRASH", y_i)
                raise
            X_i, y_i = X_train_pca[mask], y_i[mask]

            X_train_split, X_val, y_train_i, y_val_i = train_test_split(
                X_i, y_i,
                test_size=config.get("val_fraction", 0.1),
                random_state=42
            )

            prop_models = []

            for k in range(n_ensemble):             

                model = CatBoostRegressor(**catboost_params, random_seed=round_idx * 100 + k)
                model.fit(
                    X_train_split, y_train_i,
                    eval_set=(X_val, y_val_i),
                    early_stopping_rounds=config.get("early_stopping_rounds", 50),
                    verbose=False
                )

                prop_models.append(model)

            models[prop] = prop_models

            # --------------------------------------------------
            # Save round models (checkpoint)
            # --------------------------------------------------
            if config.get("save_round_models", True):
                model_dir = data_dir / "models" / f"round_{round_idx}"
                model_dir.mkdir(parents=True, exist_ok=True)

                for prop, model_list in models.items():
                    for m_idx, m in enumerate(model_list):
                        m.save_model(
                            str(model_dir / f"{prop.replace('/', '_')}_model_{m_idx}.cbm")
                        )

                # Save PCA
                
                joblib.dump(pca, model_dir / "pca.joblib")
                logger.debug(f"Saved round {round_idx} models to {model_dir}")

        # ------------------------
        # Bootstrapped evaluation
        # ------------------------

        eval_dir = plots_dir / "evaluation" / f"round{round_idx}"
        df_eval = context["df_seeds"][["sequence", "ancestor_id", "source"] + multi_condition_props].copy()
        df_eval[multi_condition_props] = df_eval[multi_condition_props].apply(pd.to_numeric, errors="coerce")

        eval_result = evaluate_with_bootstrap(
            models=models,
            pca=pca,
            embed_fn=esm_embed_cached,
            df_eval=df_eval,
            properties=multi_condition_props,
            n_boot=config.get("eval_bootstrap", 500)
        )

        save_evaluation(
            eval_result,
            eval_dir,
            score_config=config.get("scoring", None),
            train_stats=df_labeled[multi_condition_props].agg(["mean", "std"])
        )

        logger.info(
            f"[Eval R{round_idx}] "
            + ", ".join(
                f"{p}: r={eval_result.metrics[p]['pearson'].mean:.2f} "
                f"[{eval_result.metrics[p]['pearson'].low:.2f}, "
                f"{eval_result.metrics[p]['pearson'].high:.2f}]"
                for p in eval_result.metrics
            )
        )

        # ------------------------
        # Candidate generation
        # ------------------------
        logger.info(f"Generating candidate mutants for {len(current_seeds)} seeds...")
        seed_id_to_ancestor = dict(zip(current_seed_ids, current_ancestor_ids))
        candidate_seqs, candidate_meta = [], []

        for sid, seq in zip(current_seed_ids, current_seeds):
            muts = generate_mutants(seq)
            candidate_seqs.extend(muts)
            candidate_meta.extend([
                {"seed_id": sid, "ancestor_id": seed_id_to_ancestor[sid], "parent_round": round_idx} for _ in muts
            ])

        df_candidates = pd.DataFrame(candidate_meta)
        df_candidates["sequence"] = candidate_seqs
        df_candidates["round"] = round_idx
        ancestor_to_dotm = dict(
            zip(context["df_seeds"]["ancestor_id"], context["df_seeds"]["vhh_num"])
        )

        df_candidates["vhh_num"] = df_candidates["ancestor_id"].map(ancestor_to_dotm)

        # Ensure mutations are only in CDRs
        seed_id_to_parent = dict(zip(current_seed_ids, current_seeds))
        candidate_array = np.array([list(seq) for seq in candidate_seqs])
        parent_array = np.array([list(seed_id_to_parent[sid]) for sid in df_candidates["seed_id"]])
        mut_mask = candidate_array != parent_array
        allowed_mask = np.zeros(candidate_array.shape[1], dtype=bool)
        allowed_mask[mutable_positions] = True
        illegal_mut_mask = mut_mask & (~allowed_mask)
        illegal_positions = np.where(illegal_mut_mask)
        if illegal_positions[0].size > 0:
            idx = illegal_positions[0][0]
            pos = illegal_positions[1][0]
            sid = df_candidates.iloc[idx]["seed_id"]
            seq = candidate_seqs[idx]
            raise RuntimeError(
                f"Illegal mutation outside CDR at position {pos} "
                f"(ancestor_id={seed_id_to_ancestor[sid]}, seq={seq})"
            )

        # Mutation counts
        seed_enc = {sid: np.array([aa_to_int[aa] for aa in s]) for sid, s in zip(current_seed_ids, current_seeds)}
        mut_counts = []
        for sid, df_g in df_candidates.groupby("seed_id"):
            seqs_arr = np.array([[aa_to_int[aa] for aa in s] for s in df_g["sequence"]])
            counts = np.sum(seqs_arr != seed_enc[sid], axis=1)
            mut_counts.extend(counts)
        df_candidates["mut_count"] = mut_counts

        # Embed candidates and predict properties
        X_emb = esm_embed_cached(df_candidates["sequence"].tolist())
        X_cand_pca = pca.transform(X_emb)

        pred_mean, pred_std = {}, {}
        for prop, model_list in models.items():
            all_preds = np.column_stack([m.predict(X_cand_pca) for m in model_list])
            pred_mean[prop] = all_preds.mean(axis=1)
            pred_std[prop] = all_preds.std(axis=1)

        for prop in multi_condition_props:
            df_candidates[prop] = pred_mean[prop]
            df_candidates[f"{prop}_pred_std"] = pred_std[prop]
            df_candidates[f"{prop}_pred_lci"] = pred_mean[prop] - 1.96 * pred_std[prop]
            df_candidates[f"{prop}_pred_uci"] = pred_mean[prop] + 1.96 * pred_std[prop]

        # Pareto and scoring

        logger.info("Selecting top candidates based on score and Pareto fronts...")

        # ==========================================================
        # Multi-strategy scoring
        # ==========================================================

        df_candidates["pareto_flag"] = compute_multi_condition_pareto(
            df_candidates, multi_condition_props
        )

        scoring_cfg = config.get("scoring", {})
        method = scoring_cfg.get("method", "simple")
        mutation_penalty = scoring_cfg.get("mutation_penalty", 0.1)
        opt_direction = scoring_cfg.get("optimisation_direction", {})
        prop_weights = scoring_cfg.get("property_weights", {})
        use_novelty = scoring_cfg.get("use_novelty_bonus", False)
        novelty_weight = scoring_cfg.get("novelty_weight", 0.05)

        # ------------------------------------------------------------------
        # SIMPLE (original behaviour)
        # ------------------------------------------------------------------
        if method == "simple":
            df_candidates["score"] = (
                df_candidates["Tm1 (°C)"]
                - mutation_penalty * df_candidates["mut_count"]
            )

        # ------------------------------------------------------------------
        # Z-SCORED MULTI-PROPERTY UTILITY (recommended)
        # ------------------------------------------------------------------
        elif method == "zscore_weighted":

            train_stats = df_labeled[multi_condition_props].agg(["mean", "std"])

            utility = np.zeros(len(df_candidates))

            for prop in multi_condition_props:
                if prop not in df_candidates:
                    continue

                mu = train_stats.loc["mean", prop]
                sigma = train_stats.loc["std", prop] + 1e-8

                z = (df_candidates[prop] - mu) / sigma

                direction = opt_direction.get(prop, 1)
                weight = prop_weights.get(prop, 1.0)

                utility += weight * direction * z

            df_candidates["score"] = utility - mutation_penalty * df_candidates["mut_count"]

        # ------------------------------------------------------------------
        # PARETO-RANK SCORING
        # ------------------------------------------------------------------
        elif method == "pareto_rank":

            # Percentile rank across properties
            ranked = []

            for prop in multi_condition_props:
                direction = opt_direction.get(prop, 1)
                r = df_candidates[prop].rank(pct=True)

                if direction == -1:
                    r = 1 - r

                ranked.append(r)

            mean_rank = np.vstack(ranked).mean(axis=0)

            df_candidates["score"] = (
                2.0 * df_candidates["pareto_flag"].astype(float)
                + mean_rank
                - mutation_penalty * df_candidates["mut_count"]
            )

        else:
            raise ValueError(f"Unknown scoring method: {method}")

        # ------------------------------------------------------------------
        # Optional novelty bonus
        # ------------------------------------------------------------------
        if use_novelty:
            centroid = X_train_pca.mean(axis=0)
            novelty = np.linalg.norm(X_cand_pca - centroid, axis=1)
            df_candidates["score"] += novelty_weight * novelty

        # Ancestor-aware selection
        selected = [df_a.sort_values("score", ascending=False).head(min_per_ancestor)
                    for _, df_a in df_candidates.groupby("ancestor_id")]
        df_selected = pd.concat(selected)
        remaining = top_k - len(df_selected)
        if remaining > 0:
            filler = df_candidates.drop(df_selected.index).sort_values("score", ascending=False).head(remaining)
            df_selected = pd.concat([df_selected, filler])
        df_selected = df_selected.sort_values("score", ascending=False).head(top_k)
        df_selected["round"] = round_idx

        # Update seeds for next round
        current_seeds = df_selected["sequence"].tolist()
        current_seed_ids = list(range(len(current_seeds)))
        current_ancestor_ids = df_selected["ancestor_id"].tolist()

        # Simulated measurement noise
        df_measured = df_selected.copy()

        for prop in multi_condition_props:
            # Ensure numeric before adding noise
            df_measured[prop] = pd.to_numeric(df_measured[prop], errors="coerce").astype(float)
            # Add Gaussian noise
            df_measured[prop] += np.random.normal(0, noise_scale, size=len(df_measured))

        # Concatenate into df_labeled
        df_labeled_numeric = pd.concat(
            [df_labeled, df_measured[["sequence", "ancestor_id", "vhh_num"] + multi_condition_props]],
            ignore_index=True
        )

        # Final enforcement: ensure all multi-condition columns are float
        df_labeled_numeric[multi_condition_props] = df_labeled_numeric[multi_condition_props].apply(pd.to_numeric, errors="coerce")

        # Update y_train
        y_train = df_labeled_numeric[multi_condition_props].values
        X_train_full = esm_embed_cached(df_labeled_numeric["sequence"].tolist())

        # Ensure all multi-condition columns are numeric
        for prop in multi_condition_props:
            df_labeled[prop] = pd.to_numeric(df_labeled[prop], errors="coerce")

        # Save CSVs
        df_candidates.to_csv(data_dir / f"round{round_idx}_candidates.csv", index=False)
        df_selected.to_csv(data_dir / f"round{round_idx}_selected_batch.csv", index=False)

        df_selected["round"] = round_idx
        df_selected.to_csv(data_dir / f"round{round_idx}_selected_batch.csv", index=False)

        # UMAP projection
        
        n_umap_sample = min(1000, len(df_candidates))
        sample_idx = np.random.choice(len(df_candidates), n_umap_sample, replace=False)
        reducer = umap.UMAP(n_components=2,
                            n_neighbors=config.get("umap_n_neighbors", 15),
                            min_dist=config.get("umap_min_dist", 0.1),
                            n_epochs=config.get("umap_n_epochs", 50),
                            metric="euclidean",
                            n_jobs=-1)
        X_umap_sample = reducer.fit_transform(X_emb[sample_idx])

        # Nearest-neighbor mapping for all candidates
        nn = NearestNeighbors(n_neighbors=1).fit(X_emb[sample_idx])
        _, idx = nn.kneighbors(X_emb)
        df_candidates["umap1"] = X_umap_sample[idx[:, 0], 0]
        df_candidates["umap2"] = X_umap_sample[idx[:, 0], 1]

        # Color mapping for ancestors
        surviving_ancestors = sorted(df_candidates["ancestor_id"].unique())
        tab_palettes = [cm.tab20.colors, cm.tab20b.colors, cm.tab20c.colors]
        combined_colours = np.vstack(tab_palettes)
        if round_idx == 1:
            ancestor_colour_map = {anc: combined_colours[i % len(combined_colours)]
                                   for i, anc in enumerate(surviving_ancestors)}

        # UMAP plot

        logger.info("Generating UMAP and radar plots...")

        plt.figure(figsize=(9, 7))
        sns.scatterplot(data=df_candidates, x="umap1", y="umap2",
                        color="lightgray", s=25, alpha=0.35, legend=False)
        df_pareto = df_candidates[df_candidates["pareto_flag"]]
        sns.scatterplot(data=df_pareto, x="umap1", y="umap2",
                        hue="ancestor_id", palette=ancestor_colour_map,
                        s=90, alpha=0.95, legend="full")
        plt.title(f"Round {round_idx} UMAP – Pareto Regions by Ancestor")
        plt.tight_layout()
        plt.savefig(plots_dir / f"round{round_idx}_umap_pareto_ancestors.png", dpi=150)
        plt.close()

        # Radar plot
        plot_round_radar(df_candidates, round_idx, plots_dir / f"round{round_idx}_radar.png", multi_condition_props)

        round_elapsed = time.perf_counter() - round_start_time

        logger.info(f"Round {round_idx} stats: median Tm={df_candidates['Tm1 (°C)'].median():.2f}, "
            f"mean score={df_candidates['score'].mean():.2f}, "
            f"Pareto count={int(df_candidates['pareto_flag'].sum())}")
        
        logger.info(
            f"Round {round_idx} timing: "
            f"{round_elapsed:.1f}s "
            f"({round_elapsed/60:.2f} min)"
        )
        
        logger.info(f"Round {round_idx} complete.\n")        
        
        round_stat = {
            "round": round_idx,
            "mean_score": df_candidates["score"].mean(),
            "pareto_count": int(df_candidates["pareto_flag"].sum())
        }

        for prop in multi_condition_props:

            if prop not in df_candidates.columns:
                continue

            # Median experimental
            round_stat[f"median_{prop}"] = df_candidates[prop].median()
            round_stat[f"mean_{prop}"]   = df_candidates[prop].mean()

            # --- Predicted uncertainty aggregation ---
            lci_col = f"{prop}_pred_lci"
            uci_col = f"{prop}_pred_uci"

            if lci_col in df_candidates.columns:
                round_stat[f"median_{prop}_lci"] = df_candidates[lci_col].median()
                round_stat[f"median_{prop}_uci"] = df_candidates[uci_col].median()

            # --- Pareto-only summaries ---
            df_p = df_candidates[df_candidates["pareto_flag"]]

            if not df_p.empty:
                round_stat[f"pareto_median_{prop}"] = df_p[prop].median()
                round_stat[f"pareto_mean_{prop}"]   = df_p[prop].mean()

        context["round_stats"].append(round_stat)

    context.update({
        "df_labeled": df_labeled_numeric,
        "current_seeds": current_seeds,
        "current_seed_ids": current_seed_ids,
        "current_ancestor_ids": current_ancestor_ids,
        "esm_cache": esm_cache,
        "pca": pca
    })

    # Summary plots and CSVs

    logger.info("Generating summary plots and evaluations...")

    cdr_regions = context["cdr_regions"]

    generate_summary_plots(
        df_all=pd.concat(
            [pd.read_csv(data_dir / f"round{r}_candidates.csv")
            for r in range(1, n_rounds + 1)],
            ignore_index=True
        ),
        centroid_history=centroid_history,
        ancestor_colour_map=ancestor_colour_map,
        context=context,
        data_dir=data_dir,
        plots_dir=plots_dir,
        ancestor_to_seq=dict(zip(
            context["df_seeds"]["ancestor_id"],
            context["df_seeds"]["sequence"]
        )),
        aa_to_int={aa: i for i, aa in enumerate("ACDEFGHIKLMNPQRSTVWY")},
        cdr_regions=cdr_regions,
        framework_positions=context["framework_positions"],
        max_variants=config.get("max_variants_per_heatmap", 150),
    )

    # ============================================================
    # Final evaluation + executive summary
    # ============================================================

    logger.info("Evaluation: running round-wise evaluation with bootstrapping...")

    df_all_rounds = pd.concat(
        [pd.read_csv(data_dir / f"round{r}_candidates.csv") for r in range(1, n_rounds + 1)],
        ignore_index=True
    )

    # ------------------------------------------------
    # Precompute embeddings + PCA for all sequences
    # ------------------------------------------------
    logger.info("Precomputing embeddings for round-wise evaluation...")

    all_sequences = df_all_rounds["sequence"].tolist()
    X_all_emb = esm_embed_cached(all_sequences)
    X_all_pca = pca.transform(X_all_emb)

    eval_dir = plots_dir / "evaluation"

    df_metrics = run_roundwise_evaluation(
        models=models,
        pca=pca,
        embed_fn=esm_embed_cached,
        X_all_pca = X_all_pca,
        df_all=df_all_rounds,
        properties=multi_condition_props,
        out_dir=eval_dir,
        n_boot=500,
        min_n=20,
    )

    plot_learning_curves_per_property(
        df_metrics,
        out_path=eval_dir / "executive_learning_curves.png"
    )

    logger.info("All plots and evaluations complete")

    logger.info("Active learning rounds complete")

    # ============================================================
    # Train final production model on ALL labeled data
    # ============================================================

    if config.get("train_final_model", True):

        logger.info("Training final production models on full dataset...")

        # Re-embed all labeled sequences
        df_full = context["df_labeled"].copy()
        X_full_emb = esm_embed_cached(df_full["sequence"].tolist())
        X_full_pca = pca.transform(X_full_emb)

        final_models = {}
        n_ensemble = config.get("n_model_ensemble", 5)

        for i, prop in enumerate(multi_condition_props):

            y = pd.to_numeric(df_full[prop], errors="coerce").astype(float).values
            mask = np.isfinite(y)
            X_i, y_i = X_full_pca[mask], y[mask]

            prop_models = []

            for k in range(n_ensemble):
                model = CatBoostRegressor(
                    **catboost_params,
                    random_seed=9999 + k
                )
                model.fit(X_i, y_i, verbose=False)
                prop_models.append(model)

            final_models[prop] = prop_models

        # Save final production bundle
        prod_dir = data_dir / "models" / "production"
        prod_dir.mkdir(parents=True, exist_ok=True)

        for prop, model_list in final_models.items():
            for m_idx, m in enumerate(model_list):
                m.save_model(
                    str(prod_dir / f"{prop.replace('/', '_')}_model_{m_idx}.cbm")
                )

        joblib.dump(pca, prod_dir / "pca.joblib")

        import json
        with open(prod_dir / "metadata.json", "w") as f:
            json.dump({
                "properties": multi_condition_props,
                "n_ensemble": n_ensemble,
                "catboost_params": catboost_params,
            }, f, indent=2)

        logger.info(f"Saved final production model to {prod_dir}")
        
    return context