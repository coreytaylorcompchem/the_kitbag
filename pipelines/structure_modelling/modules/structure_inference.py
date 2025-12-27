import json
import csv

from pathlib import Path

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
