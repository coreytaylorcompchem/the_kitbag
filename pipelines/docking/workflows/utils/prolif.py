import re

from pathlib import Path
import numpy as np
import pandas as pd
from collections import OrderedDict
from tqdm import tqdm

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from rdkit import Chem

import prolif as plf
from prolif import sdf_supplier

from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

def natural_key(s):
    return [
        int(text) if text.isdigit() else text.lower()
        for text in re.split(r"(\d+)", str(s))
    ]

def select_best_poses(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)

    # Defensive check: only use docks that succeed.
    if "dock_status" in df.columns:
        df = df[df["dock_status"] == "success"]

    df = df.dropna(subset=["docking_score"])
    df = df[df["docking_score"] < 0]

    if df.empty:
        raise ValueError("No successful docking poses found after filtering.")

    idx = df.groupby("name")["docking_score"].idxmin()
    best_df = df.loc[idx].reset_index(drop=True)

    return best_df



def build_prolif_dataframe(
    best_pose_df: pd.DataFrame,
    protein_pdb: Path,
) -> pd.DataFrame:

    protein_rdkit = Chem.MolFromPDBFile(
        str(protein_pdb),
        removeHs=False
    )
    protein = plf.Molecule.from_rdkit(protein_rdkit)

    dfs = []
    ligand_names = []

    for _, row in tqdm(
            best_pose_df.iterrows(),
            total=len(best_pose_df),
            desc="Calculating ProLIF fingerprints",
        ):
        sdf_path = row.get("sdf_path")
        pose_idx = row.get("pose_idx")

        if not sdf_path or not Path(sdf_path).exists():
            continue
        if pose_idx is None:
            continue

        ligands = list(sdf_supplier(str(sdf_path)))
        if pose_idx >= len(ligands):
            continue

        mol = ligands[int(pose_idx)]
        if mol is None:
            logger.warning(
                f"Skipping ligand {row['name']} – RDKit failed to parse SDF pose"
            )
            continue

        ligand = plf.Molecule.from_rdkit(mol)

        try:
            fp = plf.Fingerprint()
            fp.run_from_iterable([ligand], protein, progress=False)
        except Exception as e:
            logger.warning(
                f"ProLIF failed for ligand {row['name']}: {e}"
            )
            continue

        df_fp_single = fp.to_dataframe()
        dfs.append(df_fp_single)
        ligand_names.append(row["name"])

    if not dfs:
        raise ValueError("No valid ProLIF fingerprints generated.")

    # Build DataFrame directly avoiding Prolif's crap helpers

    df_fp = pd.concat(dfs, axis=0)
    df_fp.index = ligand_names
    df_fp.index.name = "ligand"
    df_fp = (
        df_fp
        .astype("boolean")
        .fillna(False)
    )

    return df_fp

def plot_prolif_barcode(
    fp_df,
    output_dir: Path,
    interaction_colors: dict,
    min_frequency: float = 0.1,
):
    mat = fp_df.astype(int)

    interaction_level = mat.columns.nlevels - 1
    residue_level = interaction_level - 1

    # Filter rare interactions
    interaction_freq = (
        mat.T
        .groupby(level=interaction_level)
        .mean()
        .mean(axis=1)
    )

    keep_interactions = interaction_freq[
        interaction_freq >= min_frequency
    ].index

    mat = mat.loc[
        :,
        mat.columns.get_level_values(interaction_level)
                   .isin(keep_interactions)
    ]

    # Natural sort ligand IDs
    sorted_ligands = sorted(mat.index, key=natural_key)
    mat = mat.loc[sorted_ligands]

    mat = mat.sort_index(axis=1)
    mat_plot = mat.T  # residues × ligands

    # Build residue groups
    residues = mat_plot.index.get_level_values(residue_level)
    unique_residues = list(OrderedDict.fromkeys(residues))

    y_ticks = []
    y_labels = []
    row_map = []

    y = 0
    for res in unique_residues:
        idxs = [i for i, r in enumerate(residues) if r == res]
        start = y
        for i in idxs:
            row_map.append((y, i))
            y += 1
        end = y - 1

        y_ticks.append((start + end) / 2)
        y_labels.append(res)

    # Adaptive sizing for plot
    n_ligands = mat_plot.shape[1]

    # Thinner columns as ligand count grows
    col_width = min(1.0, max(0.2, 40 / n_ligands))

    fig_w = min(30, max(8, n_ligands * 0.25))
    fig_h = max(6, len(row_map) * 0.25)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # Background
    ax.imshow(
        np.zeros((len(row_map), n_ligands)),
        cmap="Greys",
        vmin=0,
        vmax=1,
        aspect="auto"
    )

    # Draw interactions
    for y_pos, row_idx in row_map:
        _, residue, interaction = mat_plot.index[row_idx]
        color = interaction_colors.get(interaction, "#000000")

        for x in range(n_ligands):
            if mat_plot.iloc[row_idx, x]:
                ax.add_patch(
                    Rectangle(
                        (x - col_width / 2, y_pos - 0.5),
                        col_width,
                        1,
                        facecolor=color,
                        edgecolor="none"
                    )
                )

    # Axes formatting
    ax.set_xticks(np.arange(n_ligands))
    ax.set_xticklabels(mat_plot.columns, rotation=45, ha="right")

    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels)

    ax.set_xlabel("Ligand ID")
    ax.set_ylabel("Protein residue")
    ax.set_title("Barcode (P-L interactions)")

    # Grid
    ax.set_xticks(np.arange(-0.5, n_ligands, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(row_map), 1), minor=True)
    ax.grid(which="minor", color="lightgrey", linestyle="-", linewidth=0.5)

    # Legend
    handles = [
        Rectangle((0, 0), 1, 1, color=c)
        for c in interaction_colors.values()
    ]
    labels = list(interaction_colors.keys())

    ax.legend(
        handles,
        labels,
        title="Interaction type",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        borderaxespad=0,
    )

    plt.tight_layout()
    out = output_dir / "prolif_barcode_best_poses.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
