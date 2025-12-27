import json
import csv
import shutil
import statistics
import subprocess
import numpy as np

from pathlib import Path

from scipy.cluster.hierarchy import linkage, fcluster
from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain, Residue, Atom

from pipeline.task_registry import register_task

from backends.structure_inference import StructureInferenceBackend
from backends.chai import ChaiBackend
from backends.boltz import BoltzBackend
from backends.openfold import OpenFoldBackend

from modules.utils.ranking import compute_scores

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=True, simple_format=True)

@register_task(
    "predict_structures",
    category="Structure modelling",
    description="Run structure prediction tools (Chai, Boltz, etc.) in parallel."
)
def predict_structures(backend, config, **kwargs):
    sequence = config["protein"]["sequence"]
    output_dir = config["output_dir"]

    tools = {
        "chai": ChaiBackend(sequence, output_dir, device=0),
        "boltz": BoltzBackend(sequence, output_dir, device=1),
        "openfold": OpenFoldBackend(sequence, output_dir, device=2),
    }

    inference = StructureInferenceBackend(
        sequence=sequence,
        tools=tools,
        output_dir=output_dir,
    )

    results = inference.run_all()

    backend.cache = {"predictions": results}
    return backend.cache
@register_task(
    "rank_predictions",
    category="Structure modelling",
    description="Rank AI-predicted structures by confidence metrics.",
)
def rank_predictions(backend, config, **kwargs):
    ranking_cfg = config["structure_prediction"].get("ranking", {})
    if not ranking_cfg.get("enabled", True):
        logger.info("Ranking disabled; skipping.")
        return {}

    output_dir = Path(config["output_dir"])
    ranking_dir = output_dir / "ranking"
    ranking_dir.mkdir(exist_ok=True)

    # --- load inference results ---
    predictions = backend.cache.get("predictions")
    if not predictions:
        raise ValueError("No predictions found in backend cache")

    metric_cfg = ranking_cfg["metrics"]
    aggregation = ranking_cfg.get("aggregation", "weighted_sum")
    normalize = ranking_cfg.get("normalize", True)

    # --- compute scores ---
    scores = compute_scores(
        predictions,
        metric_cfg,
        aggregation,
        normalize,
    )

    for entry, score in zip(predictions, scores):
        entry["score"] = score

    # --- sort ---
    predictions.sort(key=lambda x: x["score"], reverse=True)

    # --- write JSON ---
    json_path = ranking_dir / "ranking.json"
    with open(json_path, "w") as f:
        json.dump(predictions, f, indent=2)

    # --- write CSV ---
    csv_path = ranking_dir / "ranking.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "rank",
            "tool",
            "run_id",
            "seed",
            "device",
            "score",
            "plddt",
            "ptm",
            "iptm",
            "structure_path",
        ])

        for i, p in enumerate(predictions, 1):
            m = p["metrics"]
            writer.writerow([
                i,
                p["tool"],
                p["run_id"],
                p["seed"],
                p["device"],
                round(p["score"], 4),
                m.get("plddt"),
                m.get("ptm"),
                m.get("iptm"),
                p["structure_path"],
            ])

    logger.info(f"Ranking written to {ranking_dir}")

    # --- select best structures ---
    selected = {}

    if ranking_cfg.get("select", {}).get("best_overall", True):
        best = predictions[0]
        selected["best_overall"] = best
        logger.info(
            f"Best overall: {best['tool']} run {best['run_id']} "
            f"(score={best['score']:.3f})"
        )

    if ranking_cfg.get("select", {}).get("best_per_tool", True):
        best_per_tool = {}
        for p in predictions:
            tool = p["tool"]
            if tool not in best_per_tool:
                best_per_tool[tool] = p
        selected["best_per_tool"] = best_per_tool

    # --- cache results ---
    backend.cache["ranking"] = {
        "ranking_json": str(json_path),
        "ranking_csv": str(csv_path),
        "selected": selected,
    }

    return backend.cache["ranking"]

@register_task(
    "rank_models",
    category="Structure modelling",
    description="Rank predicted structures from multiple backends/tools using confidence scores.",
)
def rank_models(backend_tools, config, **kwargs):
    """
    Rank models across multiple backends and seeds.

    backend_tools: list of backend instances (OpenFold, Boltz, Chai, etc.)
    config: configuration dict
    """
    all_runs = []

    # --- collect all predictions from each backend ---
    for backend in backend_tools:
        cache_ranking = backend.cache.get("runs", [])
        for run in cache_ranking:
            score = run.get("metrics", {}).get("plddt", 0)  # default ranking metric
            all_runs.append({
                "tool": backend.name,
                "run_id": run.get("run_id", 0),
                "seed": run.get("seed", 0),
                "structure_path": run["structure_path"],
                "metrics": run.get("metrics", {}),
                "score": score,
            })

    if not all_runs:
        raise ValueError("No predictions found for ranking.")

    # --- sort by score descending ---
    ranked = sorted(all_runs, key=lambda r: r["score"], reverse=True)

    # --- save ranking to backend cache ---
    ranking_dir = Path(config["output_dir"]) / "ranking"
    ranking_dir.mkdir(exist_ok=True)
    ranking_json = ranking_dir / "ranking.json"

    with open(ranking_json, "w") as f:
        json.dump(ranked, f, indent=2)

    # --- store in cache for downstream consensus ---
    for backend in backend_tools:
        backend.cache["ranking"] = {
            "ranking_json": str(ranking_json),
            "results": ranked
        }

    logger.info(f"Ranked {len(ranked)} models across {len(backend_tools)} backends.")
    return backend_tools[0].cache["ranking"]

def compute_tmscore_matrix(pdb_files):
    """
    Compute pairwise TM-score matrix for a list of PDB files using TM-align.
    Returns a symmetric numpy array (1-TM-score as distance).
    """
    n = len(pdb_files)
    dist_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            cmd = ["TMalign", str(pdb_files[i]), str(pdb_files[j])]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            # Parse TM-score from TM-align output
            tm_score = None
            for line in result.stdout.splitlines():
                if line.startswith("TM-score="):
                    tm_score = float(line.split()[1])
                    break
            if tm_score is None:
                raise RuntimeError(f"TM-score not found for {pdb_files[i]} vs {pdb_files[j]}")
            dist = 1 - tm_score
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist

    return dist_matrix

def cluster_structures(pdb_files, cutoff=0.3):
    """
    Cluster PDB files based on pairwise TM-score distances.
    Returns cluster labels (1-based) and indices of cluster centroids.
    """
    dist_matrix = compute_tmscore_matrix(pdb_files)
    # Condensed distance matrix for linkage
    from scipy.spatial.distance import squareform
    condensed = squareform(dist_matrix)
    Z = linkage(condensed, method='average')
    labels = fcluster(Z, t=cutoff, criterion='distance')

    # Determine cluster centroids (closest to cluster mean)
    centroids = []
    for cluster_id in set(labels):
        indices = np.where(labels == cluster_id)[0]
        if len(indices) == 1:
            centroids.append(indices[0])
        else:
            cluster_distances = dist_matrix[np.ix_(indices, indices)]
            mean_dist = cluster_distances.mean(axis=1)
            centroids.append(indices[np.argmin(mean_dist)])
    return labels, centroids

@register_task(
    "generate_consensus",
    category="Structure modelling",
    description="Generate consensus / ensemble from ranked predictions with optional clustering.",
)
def generate_consensus(backend_tools, config, **kwargs):
    """
    Generate consensus structure from multiple tools and seeds.

    backend_tools: list of backend instances (OpenFold, Boltz, Chai, etc.)
    config: configuration dict
    """
    def average_structures(pdb_files, out_path):
        """Average atomic coordinates across a list of PDB files."""
        parser = PDBParser(QUIET=True)
        structures = [parser.get_structure(f"s{i}", pdb) for i, pdb in enumerate(pdb_files)]

        ref_structure = structures[0]
        for model in ref_structure:
            for chain in model:
                for res in chain:
                    for atom in res:
                        atom.coord = np.zeros(3)

        n_structs = len(structures)
        for structure in structures:
            for ref_model, model in zip(ref_structure, structure):
                for ref_chain, chain in zip(ref_model, model):
                    for ref_res, res in zip(ref_chain, chain):
                        for ref_atom in ref_res:
                            atom = res[ref_atom.get_name()]
                            ref_atom.coord += atom.coord / n_structs

        io = PDBIO()
        io.set_structure(ref_structure)
        io.save(out_path)
        return out_path

    consensus_cfg = config["structure_prediction"].get("consensus", {})
    if not consensus_cfg.get("enabled", False):
        logger.info("Consensus generation disabled; skipping.")
        return {}

    output_dir = Path(config["output_dir"])
    consensus_dir = output_dir / "consensus"
    ensemble_dir = consensus_dir / "ensemble"
    consensus_dir.mkdir(exist_ok=True)
    ensemble_dir.mkdir(exist_ok=True)

    # --- flatten all runs across backends ---
    all_runs = []
    for backend in backend_tools:
        ranking = backend.cache.get("ranking", {}).get("results", [])
        all_runs.extend(ranking)

    if not all_runs:
        raise ValueError("No predictions available for consensus.")

    # --- select top models ---
    selection_cfg = consensus_cfg.get("selection", {})
    mode = selection_cfg.get("mode", "top_k")
    selected = []

    if mode == "top_k":
        k = selection_cfg.get("k", 1)
        selected = sorted(all_runs, key=lambda r: r.get("score", 0), reverse=True)[:k]
    elif mode == "score_threshold":
        threshold = selection_cfg["score_threshold"]
        selected = [r for r in all_runs if r.get("score", 0) >= threshold]
    elif mode == "per_tool":
        seen = set()
        for r in all_runs:
            if r["tool"] not in seen:
                selected.append(r)
                seen.add(r["tool"])
    else:
        raise ValueError(f"Unknown consensus selection mode: {mode}")

    if not selected:
        raise ValueError("No models selected for consensus.")

    logger.info(f"Selected {len(selected)} models for consensus")

    # --- copy PDBs into ensemble dir ---
    pdb_files = []
    for model in selected:
        run_id = model.get("run_id", 0)
        seed = model.get("seed", 0)
        src = Path(model["structure_path"])
        dst = ensemble_dir / f"{model['tool']}_run{run_id}_seed{seed}.pdb"
        shutil.copy(src, dst)
        pdb_files.append(dst)
        model["ensemble_path"] = str(dst)

    # --- optional structural clustering ---
    clustering_enabled = consensus_cfg.get("clustering", {}).get("enabled", True)
    cluster_cutoff = consensus_cfg.get("clustering", {}).get("cutoff", 0.3)
    cluster_labels = None
    cluster_centroids = None

    if clustering_enabled and len(pdb_files) > 1:
        cluster_labels, centroid_indices = cluster_structures(pdb_files, cutoff=cluster_cutoff)
        cluster_centroids = [selected[i] for i in centroid_indices]
        logger.info(f"Identified {len(set(cluster_labels))} clusters")

    # --- select representative ---
    rep_cfg = consensus_cfg.get("representative", {})
    representative = None
    if rep_cfg.get("enabled", True):
        strategy = rep_cfg.get("strategy", "best_score")
        candidates = cluster_centroids if clustering_enabled and cluster_centroids else selected

        if strategy == "best_score":
            representative = max(candidates, key=lambda x: x.get("score", 0))
        elif strategy == "median_score":
            scores = [m.get("score",0) for m in candidates]
            median = statistics.median(scores)
            representative = min(candidates, key=lambda x: abs(x.get("score",0)-median))
        else:
            raise ValueError(f"Unknown representative strategy: {strategy}")

        rep_dst = consensus_dir / "representative.pdb"
        shutil.copy(Path(representative["structure_path"]), rep_dst)
        logger.info(f"Representative model: {representative['tool']} run {representative['run_id']} seed {representative.get('seed',0)} (score={representative.get('score',0):.3f})")

    # --- optional averaging of centroids ---
    average_enabled = consensus_cfg.get("average_centroids", {}).get("enabled", True)
    average_structure = None
    if average_enabled and cluster_centroids and len(cluster_centroids) > 1:
        centroid_files = [Path(c["structure_path"]) for c in cluster_centroids]
        avg_path = consensus_dir / "consensus_average.pdb"
        average_structures(centroid_files, avg_path)
        logger.info(f"Averaged consensus structure written to: {avg_path}")
        average_structure = str(avg_path)

    # --- save consensus metadata ---
    consensus_data = {
        "n_models": len(selected),
        "mean_score": sum(m.get("score",0) for m in selected) / len(selected),
        "max_score": max(m.get("score",0) for m in selected),
        "ensemble": selected,
        "cluster_labels": cluster_labels.tolist() if cluster_labels is not None else None,
        "cluster_centroids": [c["tool"] + f"_run{c['run_id']}_seed{c.get('seed',0)}" for c in cluster_centroids] if cluster_centroids else None,
        "representative": representative,
        "average_structure": average_structure,
    }

    with open(consensus_dir / "consensus.json", "w") as f:
        json.dump(consensus_data, f, indent=2)

    return {
        "consensus_dir": str(consensus_dir),
        "ensemble_dir": str(ensemble_dir),
        "representative": str(consensus_dir / "representative.pdb") if representative else None,
        "average_structure": average_structure,
        "metadata": consensus_data,
    }
