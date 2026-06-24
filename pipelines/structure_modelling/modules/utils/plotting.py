import matplotlib.pyplot as plt


def _collect_metric(predictions, key):
    return [
        p["metrics"].get(key)
        for p in predictions
        if p["metrics"].get(key) is not None
    ]


def plot_metric(predictions, key, out_path):
    values = _collect_metric(predictions, key)

    if not values:
        return

    plt.figure()
    plt.hist(values, bins=20)
    plt.xlabel(key)
    plt.ylabel("Count")
    plt.title(f"{key} distribution")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


def plot_score(predictions, out_path):
    values = [p.get("score") for p in predictions if p.get("score") is not None]

    if not values:
        return

    plt.figure()
    plt.hist(values, bins=20)
    plt.xlabel("score")
    plt.ylabel("Count")
    plt.title("Combined score distribution")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


def scatter_plot(predictions, x_key, y_key, out_path):
    x = []
    y = []

    for p in predictions:
        xv = p["metrics"].get(x_key)
        yv = p["metrics"].get(y_key)

        if xv is not None and yv is not None:
            x.append(xv)
            y.append(yv)

    if not x:
        return

    plt.figure()
    plt.scatter(x, y)
    plt.xlabel(x_key)
    plt.ylabel(y_key)
    plt.title(f"{x_key} vs {y_key}")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()