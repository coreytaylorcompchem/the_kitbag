from typing import List, Dict


def normalize_metric(values):
    clean = [v for v in values if v is not None]
    if not clean:
        return {i: None for i in range(len(values))}

    vmin, vmax = min(clean), max(clean)
    if vmin == vmax:
        return {i: 1.0 for i in range(len(values))}

    return {
        i: (v - vmin) / (vmax - vmin) if v is not None else None
        for i, v in enumerate(values)
    }


def compute_scores(
    entries: List[dict],
    metric_cfg: dict,
    aggregation: str,
    normalize: bool = True,
):
    metric_names = metric_cfg.keys()

    # collect raw values
    raw = {
        m: [e["metrics"].get(m) for e in entries]
        for m in metric_names
    }

    # normalize if needed
    norm = {}
    for m, values in raw.items():
        if normalize:
            norm[m] = normalize_metric(values)
        else:
            norm[m] = {i: v for i, v in enumerate(values)}

    scores = []
    for i, entry in enumerate(entries):
        s = 0.0
        count = 0

        for m, cfg in metric_cfg.items():
            value = norm[m].get(i)
            if value is None:
                continue

            weight = cfg.get("weight", 1.0)
            direction = 1 if cfg.get("higher_is_better", True) else -1

            s += direction * weight * value
            count += 1

        if aggregation == "mean" and count > 0:
            s /= count

        scores.append(s)

    return scores
