import random
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import esm
from catboost import CatBoostRegressor
from joblib import Parallel, delayed
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
import umap
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger
from modules.utils.plotting import (
    compute_multi_condition_pareto,
    plot_round_radar
)

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

# ---------------------------
# Task 1: load_seed_dataset
# ---------------------------

@register_task("load_seed_dataset", category="VHH", description="Load seed sequences and initialize context")
def load_seed_dataset(config, context):
    df = pd.read_csv(config["csv_path"])
    expected_len = config.get("expected_len", 125)
    df = df[df["sequence"].str.len() == expected_len].reset_index(drop=True)

    seeds = df["sequence"].tolist()
    seed_ids = list(range(len(df)))

    # define mutable positions
    seq_len = len(seeds[0])
    framework_positions = list(range(0, 30)) + list(range(45, 60)) + list(range(67, 95)) + list(range(115, 120))
    mutable_positions = list(set(range(seq_len)) - set(framework_positions))

    context.update({
        "df_seeds": df,
        "seeds": seeds,
        "seed_ids": seed_ids,
        "multi_condition_props": config["multi_condition_props"],
        "mutable_positions": mutable_positions,
        "current_seeds": seeds.copy(),
        "current_seed_ids": seed_ids.copy(),
        "current_ancestor_ids": df["ancestor_id"].tolist(),
        "round_stats": [],
        "centroid_history": defaultdict(list)
    })

    logger.info(f"Loaded {len(seeds)} seed sequences")
    return context

# ---------------------------
# Task 2: load_esm_model
# ---------------------------

@register_task("load_esm_model", category="VHH", description="Load pretrained ESM model")
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

# ---------------------------
# Task 3: active_learning_rounds
# ---------------------------

@register_task("active_learning_rounds", category="VHH", description="Run AL rounds with UMAP and radar plots")
def active_learning_rounds(config, context):
    from itertools import combinations
    import matplotlib.cm as cm

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
    df_labeled = context["df_seeds"][["sequence"] + multi_condition_props].copy()
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

    # --- Helpers ---
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

        for i in range(0, len(to_compute), batch_size):
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

    # --- PCA setup ---
    pca = None

    for round_idx in range(1, n_rounds+1):
        logger.info(f"=== Round {round_idx} ===")

        # --- Embed labeled sequences ---
        X_train_full = esm_embed_cached(df_labeled["sequence"].tolist())
        y_train = df_labeled[multi_condition_props].values

        if pca is None:
            n_samples, _ = X_train_full.shape
            pca_components = min(config.get("n_components", 512), n_samples)
            pca = PCA(n_components=pca_components, random_state=config.get("pca_seed", 42))
            pca.fit(X_train_full)
        X_train_pca = pca.transform(X_train_full)

        # --- Train CatBoost per property ---
        models = {}
        for i, prop in enumerate(multi_condition_props):
            y_i = y_train[:, i]
            mask = np.isfinite(y_i)
            X_i, y_i = X_train_pca[mask], y_i[mask]
            X_train, X_val, y_train_i, y_val_i = train_test_split(X_i, y_i, test_size=config.get("val_fraction", 0.1), random_state=42)
            model = CatBoostRegressor(**catboost_params, random_seed=round_idx)
            model.fit(X_train, y_train_i, eval_set=(X_val, y_val_i),
                      early_stopping_rounds=config.get("early_stopping_rounds", 50), verbose=False)
            models[prop] = model

        # --- Generate candidates ---
        seed_id_to_ancestor = dict(zip(current_seed_ids, current_ancestor_ids))
        candidate_seqs, candidate_meta = [], []
        for sid, seq in zip(current_seed_ids, current_seeds):
            muts = generate_mutants(seq)
            candidate_seqs.extend(muts)
            candidate_meta.extend([{"seed_id": sid, "ancestor_id": seed_id_to_ancestor[sid], "parent_round": round_idx} for _ in muts])
        df_candidates = pd.DataFrame(candidate_meta)
        df_candidates["sequence"] = candidate_seqs

        # --- Mutation counts ---
        seed_enc = {sid: np.array([aa_to_int[aa] for aa in s]) for sid, s in zip(current_seed_ids, current_seeds)}
        mut_counts = []
        for sid, df_g in df_candidates.groupby("seed_id"):
            seqs_arr = np.array([[aa_to_int[aa] for aa in s] for s in df_g["sequence"]])
            counts = np.sum(seqs_arr != seed_enc[sid], axis=1)
            mut_counts.extend(counts)
        df_candidates["mut_count"] = mut_counts

        # --- Embed candidates ---
        X_emb = esm_embed_cached(df_candidates["sequence"].tolist())
        X_cand_pca = pca.transform(X_emb)

        # --- Predict properties ---
        preds = np.column_stack([model.predict(X_cand_pca) for model in models.values()])
        df_preds = pd.DataFrame(preds, columns=multi_condition_props)
        df_candidates = pd.concat([df_candidates.reset_index(drop=True), df_preds], axis=1)

        # --- Pareto and scoring ---
        df_candidates["pareto_flag"] = compute_multi_condition_pareto(df_candidates, multi_condition_props)
        df_candidates["score"] = df_candidates["Tm1 (°C)"] - 0.5 * df_candidates["mut_count"]

        # --- Ancestor-aware selection ---
        selected = [df_a.sort_values("score", ascending=False).head(min_per_ancestor)
                    for _, df_a in df_candidates.groupby("ancestor_id")]
        df_selected = pd.concat(selected)
        remaining = top_k - len(df_selected)
        if remaining > 0:
            filler = df_candidates.drop(df_selected.index).sort_values("score", ascending=False).head(remaining)
            df_selected = pd.concat([df_selected, filler])
        df_selected = df_selected.sort_values("score", ascending=False).head(top_k)

        # --- Update seeds for next round ---
        current_seeds = df_selected["sequence"].tolist()
        current_seed_ids = list(range(len(current_seeds)))
        current_ancestor_ids = df_selected["ancestor_id"].tolist()

        # --- Simulated measurement noise ---
        df_measured = df_selected.copy()
        for prop in multi_condition_props:
            df_measured[prop] += np.random.normal(0, noise_scale, size=len(df_measured))
        df_labeled = pd.concat([df_labeled, df_measured[["sequence"] + multi_condition_props]], ignore_index=True)

        # --- Save CSVs ---
        df_candidates.to_csv(data_dir / f"round{round_idx}_candidates.csv", index=False)
        df_selected.to_csv(data_dir / f"round{round_idx}_selected_batch.csv", index=False)

        # --- UMAP projection ---
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

        # --- Color mapping for ancestors ---
        surviving_ancestors = sorted(df_candidates["ancestor_id"].unique())
        tab_palettes = [cm.tab20.colors, cm.tab20b.colors, cm.tab20c.colors]
        combined_colours = np.vstack(tab_palettes)
        if round_idx == 1:
            ancestor_colour_map = {anc: combined_colours[i % len(combined_colours)]
                                   for i, anc in enumerate(surviving_ancestors)}

        # --- UMAP plot ---
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

        # --- Radar plot ---
        plot_round_radar(df_candidates, round_idx, plots_dir / f"round{round_idx}_radar.png", multi_condition_props)

    # --- Save context ---
    context.update({
        "df_labeled": df_labeled,
        "current_seeds": current_seeds,
        "current_seed_ids": current_seed_ids,
        "current_ancestor_ids": current_ancestor_ids,
        "esm_cache": esm_cache,
        "pca": pca
    })
    logger.info("Active learning rounds complete")
    return context

