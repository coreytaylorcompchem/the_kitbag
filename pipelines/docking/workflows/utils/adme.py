import re

from pathlib import Path
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

def ligand_sort_key(name):
    """
    Extract integer index from ligand names like 'ligand12'.
    Falls back to name if no integer found.
    """
    m = re.search(r"(\d+)$", str(name))
    return int(m.group(1)) if m else name

def draw_threshold_lines(ax, thresholds, direction="higher_is_better"):
    """
    Draw horizontal threshold lines for ADME interpretation.
    Colors:
        good     -> green
        mediocre -> orange
        bad      -> red
    """
    if direction == "higher_is_better":
        ax.axhline(thresholds["bad"],      color="red",    linestyle="--", linewidth=1)
        ax.axhline(thresholds["mediocre"], color="orange", linestyle="--", linewidth=1)
        ax.axhline(thresholds["good"],     color="green",  linestyle="--", linewidth=1)
    else:
        ax.axhline(thresholds["good"],     color="green",  linestyle="--", linewidth=1)
        ax.axhline(thresholds["mediocre"], color="orange", linestyle="--", linewidth=1)
        ax.axhline(thresholds["bad"],      color="red",    linestyle="--", linewidth=1)


ADME_THRESHOLDS = {
    "caco2": {
        "bad": 5.0,
        "mediocre": 10.0,
        "good": 20.0
    },
    "logd": {
        "good": 2.0,    
        "mediocre": 3.0,
        "bad": 4.0
    },
    "cyp3a4": {
        "good": 5.0,
        "mediocre": 6.0,
        "bad": 7.0
    },
    "herg": {
        "good": 0.3,
        "mediocre": 0.6,
        "bad": 1.0
    }
}

DIRECTION = {
        "caco2": "higher_is_better",
        "logd": "lower_is_better",
        "cyp3a4": "lower_is_better",
        "herg": "lower_is_better",
    }


DESCRIPTORS = ["MW", "QED", "HBD", "csp3", "TPSA", "LogP"]


def plot_adme_descriptor_grids(csv_path: Path, output_dir: Path):
    df = pd.read_csv(csv_path)

    # one row per ligand
    df = (
        df.groupby("name", as_index=False)
          .first()
    )

    # enforce ligand order
    df["ligand_order"] = df["name"].apply(ligand_sort_key)
    df = df.sort_values("ligand_order").reset_index(drop=True)
    df["mol_idx"] = df["ligand_order"]

    for adme, threshold in ADME_THRESHOLDS.items():
        if adme not in df.columns:
            continue

        fig, axes = plt.subplots(3, 2, figsize=(12, 10), sharex=True)
        axes = axes.flatten()

        for ax, desc in zip(axes, DESCRIPTORS):
            if desc not in df.columns:
                ax.axis("off")
                continue

            sc = ax.scatter(
                df["mol_idx"],
                df[adme],
                c=df[desc],
                cmap="inferno_r",
                s=40,
                alpha=0.8
            )

    for adme, thresholds in ADME_THRESHOLDS.items():
        if adme not in df.columns:
            continue

        fig, axes = plt.subplots(3, 2, figsize=(12, 10), sharex=True)
        axes = axes.flatten()

        for ax, desc in zip(axes, DESCRIPTORS):
            if desc not in df.columns:
                ax.axis("off")
                continue

            sc = ax.scatter(
                df["mol_idx"],
                df[adme],
                c=df[desc],
                cmap="inferno_r",
                s=40,
                alpha=0.8
            )

            ax.set_title(f"{adme.upper()} coloured by {desc}")
            ax.set_ylabel(adme)

            # draw threshold lines
            draw_threshold_lines(
                ax,
                thresholds,
                direction=DIRECTION.get(adme, "lower_is_better")
            )

            cbar = plt.colorbar(sc, ax=ax)
            cbar.set_label(desc)

        for ax in axes:
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))

        for ax in axes[-2:]:
            ax.set_xlabel("Molecule number")

        plt.suptitle(
            f"Predicted {adme.upper()} by molecule.",
            fontsize=14
        )
        plt.tight_layout(rect=[0, 0, 1, 0.96])

        out = output_dir / f"{adme}_descriptor_grid.png"
        plt.savefig(out, dpi=300)
        plt.close()