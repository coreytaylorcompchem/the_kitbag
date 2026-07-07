import json
import csv
import shutil
import statistics
import subprocess
import numpy as np

from pathlib import Path
from collections import defaultdict

from concurrent.futures import ThreadPoolExecutor, as_completed

from scipy.cluster.hierarchy import linkage, fcluster

from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain, Residue, Atom
from Bio.PDB import MMCIFParser

from scipy.spatial.distance import squareform

parser = MMCIFParser(QUIET=True)

from pipeline.task_registry import register_task

# # from backends.structure_inference import StructureInferenceBackend
# # from backends.chai import ChaiBackend
from backends.boltz import BoltzBackend
from backends.intellifold import IntelliFoldBackend
# # from backends.openfold import OpenFoldBackend

from modules.utils.ranking import compute_scores
from modules.utils.csv_loader import load_sequences
# from modules.utils.plotting import plot_metric, plot_score, scatter_plot, plot_clusters
from modules.utils.msas import generate_local_msa_for_sequence, sequence_hash

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def has_metric(predictions, metric):
    return any(
        p["metrics"].get(metric) is not None
        for p in predictions
    )

def metric_or_default(p, metric_name, default=-1.0):
    value = p["metrics"].get(metric_name)
    return value if value is not None else default


def make_selection_tag(p):
    iptm = p["metrics"].get("iptm")
    plddt = p["metrics"].get("plddt")
    affinity_probability = p["metrics"].get("affinity_probability_binary")
    affinity_value = p["metrics"].get("affinity_pred_value")
    ranking_score = p["metrics"].get("ranking_score")

    if iptm is not None and iptm > 0:
        return f"iptm{iptm:.3f}"

    if affinity_probability is not None:
        return f"affprob{affinity_probability:.3f}"

    if affinity_value is not None:
        return f"affval{affinity_value:.3f}"

    if plddt is not None:
        return f"plddt{plddt:.3f}"

    if ranking_score is not None:
        return f"rankscore{ranking_score:.3f}"

    return "scoreNA"

def is_interacting(p):
    iptm = p["metrics"]["iptm"]
    return iptm is not None and iptm > 0.0

def get_backend(tool_name: str, config):
    tool_name = tool_name.lower()

    if tool_name == "boltz":
        return BoltzBackend(config=config)

    if tool_name in {"intellifold", "intfold"}:
        return IntelliFoldBackend(config=config)

    raise ValueError(f"Unknown structure prediction tool: {tool_name}")

@register_task(
    "generate_msas",
    category="Structure modelling",
    description="Generate local cached MSAs using MMseqs2/ColabFold search."
)
def generate_msas(backend, config, **kwargs):
    sp_cfg = config["structure_prediction"]
    msa_cfg = sp_cfg.get("msa", {})

    csv_path = sp_cfg["input_csv"]
    entries = load_sequences(csv_path)

    if not entries:
        raise ValueError("No entries found in CSV for MSA generation")

    if not msa_cfg.get("enabled", False):
        logger.info("[MSA] MSA generation disabled; loading CSV entries unchanged.")
        backend.cache["entries"] = entries
        return {
            "entries": entries,
            "msa_cache": None,
        }

    method = msa_cfg.get("method", "colabfold_local")

    if method != "colabfold_local":
        raise ValueError(
            f"Unsupported MSA method: {method}"
        )

    database_dir = Path(msa_cfg["database_dir"]).expanduser().resolve()
    cache_dir = Path(msa_cfg.get("cache_dir", "msa_cache"))

    if not cache_dir.is_absolute():
        cache_dir = Path(config["output_dir"]).parent / cache_dir

    cache_dir = cache_dir.resolve()

    command = msa_cfg.get(
        "command",
        "colabfold_search",
    )

    extra_args = msa_cfg.get("extra_args", [])
    reuse_existing = msa_cfg.get("reuse_existing", True)
    overwrite = msa_cfg.get("overwrite", False)
    n_jobs = int(msa_cfg.get("n_jobs", 1))

    if n_jobs < 1:
        raise ValueError("msa.n_jobs must be >= 1")

    if not database_dir.exists():
        raise FileNotFoundError(
            f"MSA database_dir does not exist: {database_dir}"
        )

    logger.info(
        f"[MSA] Local MSA generation enabled | method={method} "
        f"| database_dir={database_dir} "
        f"| cache_dir={cache_dir} "
        f"| n_jobs={n_jobs}"
    )

    manifest = {
        "method": method,
        "database_dir": str(database_dir),
        "cache_dir": str(cache_dir),
        "n_jobs": n_jobs,
        "entries": {},
    }

    # ------------------------------------------------------------------
    # Build unique sequence jobs.
    #
    # If the same sequence appears in multiple inputs, or same input with
    # different ligands/templates, only generate/search it once.
    # ------------------------------------------------------------------

    unique_jobs = {}
    assignments = []

    for entry in entries:
        entry_id = entry["id"]
        proteins = entry.get("proteins", {})

        entry.setdefault("msas", {})
        manifest["entries"][entry_id] = {}

        if not proteins:
            logger.info(
                f"[MSA] {entry_id}: no protein chains; skipping."
            )
            continue

        for chain_id, sequence in proteins.items():
            seq_hash = sequence_hash(sequence)
            sequence_id = f"{entry_id}_{chain_id}"

            assignments.append({
                "entry_id": entry_id,
                "chain_id": chain_id,
                "sequence_hash": seq_hash,
            })

            if seq_hash not in unique_jobs:
                unique_jobs[seq_hash] = {
                    "sequence_id": sequence_id,
                    "sequence": sequence,
                    "database_dir": database_dir,
                    "cache_dir": cache_dir,
                    "command": command,
                    "extra_args": extra_args,
                    "reuse_existing": reuse_existing,
                    "overwrite": overwrite,
                }

    logger.info(
        f"[MSA] Total chain assignments: {len(assignments)}"
    )

    logger.info(
        f"[MSA] Unique sequences requiring cache lookup/search: {len(unique_jobs)}"
    )

    # ------------------------------------------------------------------
    # Run unique MSA jobs in parallel.
    # ------------------------------------------------------------------

    msa_paths_by_hash = {}

    if n_jobs == 1:
        logger.info("[MSA] Running MSA jobs serially")

        for seq_hash, job in unique_jobs.items():
            msa_path = generate_local_msa_for_sequence(**job)
            msa_paths_by_hash[seq_hash] = str(msa_path)

    else:
        logger.info(
            f"[MSA] Running MSA jobs in parallel with n_jobs={n_jobs}"
        )

        with ThreadPoolExecutor(max_workers=n_jobs) as executor:
            future_to_hash = {}

            for seq_hash, job in unique_jobs.items():
                future = executor.submit(
                    generate_local_msa_for_sequence,
                    **job,
                )
                future_to_hash[future] = seq_hash

            completed = 0
            total = len(future_to_hash)

            for future in as_completed(future_to_hash):
                seq_hash = future_to_hash[future]

                try:
                    msa_path = future.result()
                except Exception as e:
                    logger.error(
                        f"[MSA] Failed MSA generation for sequence hash {seq_hash}: {e}"
                    )
                    raise

                msa_paths_by_hash[seq_hash] = str(msa_path)

                completed += 1

                logger.info(
                    f"[MSA] Completed {completed}/{total} MSA jobs"
                )

    # ------------------------------------------------------------------
    # Map MSA paths back onto entries.
    # ------------------------------------------------------------------

    for entry in entries:
        entry_id = entry["id"]
        proteins = entry.get("proteins", {})
        entry.setdefault("msas", {})

        for chain_id, sequence in proteins.items():
            seq_hash = sequence_hash(sequence)

            msa_path = msa_paths_by_hash.get(seq_hash)

            if msa_path is None:
                raise RuntimeError(
                    f"No MSA path found for entry={entry_id}, chain={chain_id}, "
                    f"sequence_hash={seq_hash}"
                )

            entry["msas"][chain_id] = msa_path

            manifest["entries"][entry_id][chain_id] = {
                "sequence_hash": seq_hash,
                "msa": msa_path,
            }

    manifest_path = cache_dir / "msa_manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    logger.info(
        f"[MSA] Wrote MSA manifest: {manifest_path}"
    )

    backend.cache["entries"] = entries
    backend.cache["msa_manifest"] = str(manifest_path)

    return {
        "entries": entries,
        "msa_manifest": str(manifest_path),
    }

@register_task(
    "predict_structures",
    category="Structure modelling",
    description="Run structure prediction tools (Chai, Boltz, etc.) in parallel."
)
def predict_structures(backend, config, **kwargs):

    entries = backend.cache.get("entries")

    if entries is None:
        csv_path = config["structure_prediction"]["input_csv"]
        entries = load_sequences(csv_path)

    if not entries:
        raise ValueError("No sequences found in CSV")

    # tools = {
    #     "chai": ChaiBackend(sequence, output_dir, device=0),
    #     "boltz": BoltzBackend(sequence, output_dir, device=1),
    #     "openfold": OpenFoldBackend(sequence, output_dir, device=2),
    # }

    # inference = StructureInferenceBackend(
    #     sequence=sequence,
    #     tools=tools,
    #     output_dir=output_dir,
    # )
    
    tools = config["structure_prediction"].get("tools", ["boltz"])
    devices_cfg = config["structure_prediction"].get("devices", {})

    all_results = []

    output_dir = Path(config["output_dir"])

    for tool_name in tools:
        backend_instance = get_backend(tool_name, config)

        # default all tools to GPU 0 unless specified
        device = devices_cfg.get(tool_name, 0)

        logger.info(f"[Predict] Running tool={tool_name} on device={device}")

        for i, entry in enumerate(entries):

            tool_output_dir = (
                output_dir /
                entry["id"] /
                tool_name
            )

            logger.info(
                f"[Predict] {entry['id']} | tool={tool_name} "
                f"| chains={list(entry['proteins'].keys())} "
                f"| msas={list(entry.get('msas', {}).keys())}"
            )

            result = backend_instance.run(
                run_id=i,
                device=device,
                output_dir=tool_output_dir,
                sequences=entry["proteins"],
                ligands=entry["ligands"],
                templates=entry.get("templates", []),
                msas=entry.get("msas", {}),
                row_constraints=entry.get("constraints", []),
            )

            for sample in result["results"]:
                logger.debug(
                    f"[{tool_name}] sample: {sample['structure']} "
                    f"plddt={sample.get('plddt')}"
                )

            result["input_id"] = entry["id"]
            all_results.append(result)

    logger.debug(f"All results: {all_results}")

    logger.info(
        f"{entry['id']} | chains={list(entry['proteins'].keys())} "
        f"| ligands={len(entry['ligands'])}"
    )

    backend.cache["predictions"] = all_results
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

    # load inference results
    runs = backend.cache.get("predictions", [])

    #DEBUG
    
    for run in runs:
        for sample in run["results"]:
            logger.debug(
                f"sample path={sample['structure']}"
            )


    predictions = []

    for run in runs:
        for sample in run["results"]:
            predictions.append({
                "tool": run["tool"],
                "structure_path": sample["structure"],              
                "metrics": {
                    "confidence_score": sample.get("confidence_score"),
                    "plddt": sample.get("plddt"),
                    "ptm": sample.get("ptm"),
                    "iptm": sample.get("iptm"),
                    "iplddt": sample.get("iplddt"),
                    "ranking_score": sample.get("ranking_score"),
                    "ligand_iptm": sample.get("ligand_iptm"),
                    "protein_iptm": sample.get("protein_iptm"),
                    "complex_pde": sample.get("complex_pde"),
                    "complex_ipde":sample.get("complex_ipde"),
                    "affinity_pred_value": sample.get("affinity_pred_value"),
                    "affinity_probability_binary": sample.get("affinity_probability_binary"),
                    "affinity_pred_value1": sample.get("affinity_pred_value1"),
                    "affinity_probability_binary1": sample.get("affinity_probability_binary1"),
                    "affinity_pred_value2": sample.get("affinity_pred_value2"),
                    "affinity_probability_binary2": sample.get("affinity_probability_binary2"),
                },
                "run_id": run["run_id"],
                "input_id": run.get("input_id"),
                "device": run["device"],
                "seed": None,
            })
    all_metric_names = sorted({
        metric
        for p in predictions
        for metric, value in p["metrics"].items()
        if value is not None
    })

    logger.debug(
        f"[Ranking] Metrics with at least one non-null value: {all_metric_names}"
    )
    
    if not runs:
        raise ValueError("No runs found in backend cache - predict_structures failed")

    if not predictions:
        raise ValueError("Runs found but no samples parsed - parsing failed")


    metric_cfg = ranking_cfg["metrics"]
    aggregation = ranking_cfg.get("aggregation", "weighted_sum")
    normalize = ranking_cfg.get("normalize", True)

    # compute scores
    scores = compute_scores(
        predictions,
        metric_cfg,
        aggregation,
        normalize,
    )

    for entry, score in zip(predictions, scores):
        entry["score"] = score

    # sort
    predictions.sort(key=lambda x: x["score"], reverse=True)

    # clustering in metric space
    # build feature matrix

    # Check to make sure the metrics are full of data
    
    for p in predictions:
        logger.debug(f"Metrics: {p['metrics']}")

    X = []
    valid_idx = []

    for i, p in enumerate(predictions):
        m = p["metrics"]

        plddt = m.get("plddt")
        ptm = m.get("ptm")
        iptm = m.get("iptm")
        iplddt = m.get("iplddt")

        plddt_str = f"{plddt:.3f}" if plddt is not None else "NA"

        logger.debug(
            f"{i+1}: "
            f"plddt={plddt_str} "
            f"path={p['structure_path']}"
        )

        if plddt is not None and ptm is not None and iptm is not None:
            X.append([
                plddt,
                iptm,
                iplddt if iplddt is not None else 0.0,
            ])
            valid_idx.append(i)

    if X:
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        k = min(4, len(X_scaled))  # safe default
        kmeans = KMeans(n_clusters=k, random_state=42, n_init="auto")
        labels = kmeans.fit_predict(X_scaled)

        # attach cluster labels back
        for idx, label in zip(valid_idx, labels):
            predictions[idx]["cluster"] = int(label)

        logger.debug(f"Assigned {k} clusters")
    else:
        logger.warning("No valid metric data for clustering")

    input_groups = defaultdict(list)

    for p in predictions:
        input_groups[p["run_id"]].append(p)

    input_summary = []

    for run_id, group in input_groups.items():
        iptms = [p["metrics"]["iptm"] for p in group if p["metrics"]["iptm"] is not None]
        plddts = [p["metrics"]["plddt"] for p in group if p["metrics"]["plddt"] is not None]
        ptms = [p["metrics"]["ptm"] for p in group if p["metrics"]["ptm"] is not None]

        if not iptms:
            continue

        input_summary.append({
            "run_id": run_id,
            "n_samples": len(group),
            "mean_iptm": float(np.mean(iptms)),
            "std_iptm": float(np.std(iptms)),
            "mean_plddt": float(np.mean(plddts)) if plddts else None,
            "mean_ptm": float(np.mean(ptms)) if ptms else None,
        })

    # write JSON
    json_path = ranking_dir / "ranking.json"
    with open(json_path, "w") as f:
        json.dump(predictions, f, indent=2)

    # write CSV
    csv_path = ranking_dir / "ranking.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "rank", "tool", "run_id", "seed", "device",
            "score",
            "confidence_score",
            "plddt", "ptm", "iptm", "iplddt",
            "ligand_iptm", "protein_iptm",
            "complex_pde", "complex_ipde",
            "affinity_pred_value",
            "affinity_probability_binary",
            "affinity_pred_value1",
            "affinity_probability_binary1",
            "affinity_pred_value2",
            "affinity_probability_binary2",
            "ranking_score",
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
                m.get("confidence_score"),
                m.get("plddt"),
                m.get("ptm"),
                m.get("iptm"),
                m.get("iplddt"),
                m.get("ligand_iptm"),
                m.get("protein_iptm"),
                m.get("complex_pde"),
                m.get("complex_ipde"),
                m.get("affinity_pred_value"),
                m.get("affinity_probability_binary"),
                m.get("affinity_pred_value1"),
                m.get("affinity_probability_binary1"),
                m.get("affinity_pred_value2"),
                m.get("affinity_probability_binary2"),
                m.get("ranking_score"),
                p["structure_path"],
            ])

    # plotting

    if ranking_cfg.get("plotting", {}).get("enabled", True):

        from modules.utils.plotting import (
            plot_metric,
            plot_score,
            scatter_plot,
            plot_clusters,
        )

        plot_metrics = [
            "confidence_score",
            "plddt",
            "ptm",
            "iptm",
            "iplddt",
            "ligand_iptm",
            "protein_iptm",
            "complex_pde",
            "complex_ipde",
            "affinity_pred_value",
            "affinity_probability_binary",
            "ranking_score",
        ]

        for metric in plot_metrics:
            if has_metric(predictions, metric):
                plot_metric(
                    predictions,
                    metric,
                    ranking_dir / f"{metric}_hist.png",
                )

        plot_score(predictions, ranking_dir / "score_hist.png")

        if has_metric(predictions, "plddt") and has_metric(predictions, "ptm"):
            scatter_plot(predictions, "plddt", "ptm", ranking_dir / "plddt_vs_ptm.png")

        if has_metric(predictions, "plddt") and has_metric(predictions, "iptm"):
            scatter_plot(predictions, "plddt", "iptm", ranking_dir / "plddt_vs_iptm.png")

        if has_metric(predictions, "ptm") and has_metric(predictions, "iptm"):
            scatter_plot(predictions, "ptm", "iptm", ranking_dir / "ptm_vs_iptm.png")

        if has_metric(predictions, "iptm") and has_metric(predictions, "iplddt"):
            scatter_plot(predictions, "iptm", "iplddt", ranking_dir / "iptm_vs_iplddt.png")

        if has_metric(predictions, "affinity_probability_binary") and has_metric(predictions, "iptm"):
            scatter_plot(
                predictions,
                "affinity_probability_binary",
                "iptm",
                ranking_dir / "affinity_probability_vs_iptm.png",
            )

        if has_metric(predictions, "affinity_pred_value") and has_metric(predictions, "affinity_probability_binary"):
            scatter_plot(
                predictions,
                "affinity_pred_value",
                "affinity_probability_binary",
                ranking_dir / "affinity_value_vs_probability.png",
            )

        plot_clusters(predictions, ranking_dir / "clusters.png")

    if ranking_cfg.get("plotting", {}).get("enabled", True):
        logger.info("Saved metric plots to ranking directory")
    else:
        logger.info("Plotting disabled; skipped metric plots")

    # select best structures
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
    
    #DEBUG
    logger.debug("===== RANKED STRUCTURES =====")

    for i, p in enumerate(predictions):
        plddt = p["metrics"].get("plddt")
        plddt_str = f"{plddt:.3f}" if plddt is not None else "NA"

        logger.debug(
            f"{i+1}: "
            f"plddt={plddt_str} "
            f"path={p['structure_path']}"
        )

    # select and copy structures
    select_cfg = ranking_cfg.get("select_structures", {})

    if select_cfg.get("enabled", True):

        k = select_cfg.get("k", 10)

        # split interacting vs non-interacting
        interacting = [p for p in predictions if is_interacting(p)]
        non_interacting = [p for p in predictions if not is_interacting(p)]

        if not predictions:
            logger.warning("No predictions available for selection")
            return

        # =========================
        # GLOBAL TOP-K
        # =========================
        global_dir = output_dir / "top_structures_global"
        global_dir.mkdir(exist_ok=True)

        global_paths = []

        if interacting:
            global_sorted = sorted(
                interacting,
                key=lambda x: x["metrics"]["iptm"],
                reverse=True
            )
            logger.info(f"[GLOBAL] Using iPTM ranking ({len(interacting)} structures)")
        else:
            global_sorted = sorted(
                non_interacting,
                key=lambda x: (
                    metric_or_default(x, "plddt"),
                    metric_or_default(x, "ranking_score"),
                    x.get("score", -1.0),
                ),
                reverse=True
            )
            logger.info(f"[GLOBAL] No interactions → using pLDDT")

        global_selected = global_sorted[:k]

        for i, p in enumerate(global_selected):
            src = Path(p["structure_path"])
            # iptm = p["metrics"]["iptm"]
            # plddt = p["metrics"]["plddt"]
            input_id = p.get("input_id", f"run{p['run_id']}")
            ext = src.suffix

            tag = make_selection_tag(p)

            dst = global_dir / f"{input_id}_rank{i+1}_{tag}{ext}"
            shutil.copy(src, dst)
            global_paths.append(str(dst))

        # =========================
        # PER-INPUT TOP-K
        # =========================
        per_input_dir = output_dir / "top_structures_per_input"
        per_input_dir.mkdir(exist_ok=True)

        grouped = defaultdict(list)
        for p in predictions:
            grouped[p["input_id"]].append(p)

        per_input_paths = []
        total_selected = 0

        for input_id, group in grouped.items():

            group_interacting = [p for p in group if is_interacting(p)]

            if group_interacting:
                sorted_group = sorted(
                    group_interacting,
                    key=lambda x: (
                        metric_or_default(x, "iptm"),
                        metric_or_default(x, "ranking_score"),
                        x.get("score", -1.0),
                    ),
                    reverse=True
                )
                logger.info(f"[PER-INPUT] {input_id}: using iPTM")
            else:
                sorted_group = sorted(
                    group,
                    key=lambda x: (
                        metric_or_default(x, "plddt"),
                        metric_or_default(x, "ranking_score"),
                        x.get("score", -1.0),
                    ),
                    reverse=True
                )
                logger.info(f"[PER-INPUT] {input_id}: no interaction → using pLDDT")

            top_group = sorted_group[:k]

            for j, p in enumerate(top_group):
                src = Path(p["structure_path"])
                #DEBUG
                
                logger.debug(
                        f"[COPY] rank={j+1} "
                        f"src={src.name}"
                    )

                # iptm = p["metrics"]["iptm"]
                # plddt = p["metrics"]["plddt"]
                ext = src.suffix

                tag = make_selection_tag(p)

                dst = per_input_dir / f"{input_id}_model{j+1}_{tag}{ext}"
                shutil.copy(src, dst)

                per_input_paths.append(str(dst))
                total_selected += 1

        logger.info(f"[PER-INPUT] Total selected structures: {total_selected}")

        # cache
        backend.cache["top_structures"] = {
            "global": global_paths,
            "per_input": per_input_paths,
        }

    # cache results
    backend.cache["ranking"] = {
        "ranking_json": str(json_path),
        "ranking_csv": str(csv_path),
        "selected": selected,
    }

    return backend.cache["ranking"]

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
        parser = MMCIFParser(QUIET=True)
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

    # flatten all runs across backends
    all_runs = []
    for backend in backend_tools:
        ranking = backend.cache.get("ranking", {}).get("results", [])
        all_runs.extend(ranking)

    if not all_runs:
        raise ValueError("No predictions available for consensus.")

    # select top models
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

    # copy PDBs into ensemble dir
    pdb_files = []
    for model in selected:
        run_id = model.get("run_id", 0)
        seed = model.get("seed", 0)
        src = Path(model["structure_path"])
        ext = Path(model["structure_path"]).suffix
        dst = ensemble_dir / f"{model['tool']}_run{run_id}_seed{seed}{ext}"

        shutil.copy(src, dst)
        pdb_files.append(dst)
        model["ensemble_path"] = str(dst)

    # optional structural clustering
    clustering_enabled = consensus_cfg.get("clustering", {}).get("enabled", True)
    cluster_cutoff = consensus_cfg.get("clustering", {}).get("cutoff", 0.3)
    cluster_labels = None
    cluster_centroids = None

    if clustering_enabled and len(pdb_files) > 1:
        cluster_labels, centroid_indices = cluster_structures(pdb_files, cutoff=cluster_cutoff)
        cluster_centroids = [selected[i] for i in centroid_indices]
        logger.info(f"Identified {len(set(cluster_labels))} clusters")

    # select representative
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

    # optional averaging of centroids
    average_enabled = consensus_cfg.get("average_centroids", {}).get("enabled", True)
    average_structure = None
    if average_enabled and cluster_centroids and len(cluster_centroids) > 1:
        centroid_files = [Path(c["structure_path"]) for c in cluster_centroids]
        avg_path = consensus_dir / "consensus_average.pdb"
        average_structures(centroid_files, avg_path)
        logger.info(f"Averaged consensus structure written to: {avg_path}")
        average_structure = str(avg_path)

    # save consensus metadata
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