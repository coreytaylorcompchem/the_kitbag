import re
import sys

from pathlib import Path
import numpy as np
import pandas as pd
from collections import OrderedDict
from tqdm import tqdm

import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for headless environments
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
    """
    Generate ProLIF fingerprints for the best docking poses.
    Fully defensive: skips invalid SDFs, missing poses, empty molecules.
    """

    protein_rdkit = Chem.MolFromPDBFile(str(protein_pdb), removeHs=False, sanitize=False)
    if protein_rdkit is None:
        logger.error(f"[ProLIF] Failed to parse protein PDB: {protein_pdb}")
        return None

    try:
        Chem.SanitizeMol(protein_rdkit)
    except Exception as e:
        logger.warning(f"[ProLIF] Protein sanitization failed: {e}")

    try:
        protein = plf.Molecule.from_rdkit(protein_rdkit)
    except Exception as e:
        logger.error(f"[ProLIF] Failed to create protein molecule: {e}")
        return None

    dfs = []
    ligand_names = []

    for _, row in tqdm(best_pose_df.iterrows(),
                        total=len(best_pose_df),
                        desc="Calculating ProLIF fingerprints",
                        leave=True,
                        ncols=100,
                        file=sys.stdout):

        sdf_path = row.get("sdf_path")
        pose_idx = row.get("pose_idx")
        ligand_name = row.get("name", "unknown")

        # ---- Defensive checks ----
        if not sdf_path or not Path(sdf_path).exists():
            logger.warning(f"[ProLIF] Skipping {ligand_name} – missing SDF: {sdf_path}")
            continue

        if pose_idx is None:
            logger.warning(f"[ProLIF] Skipping {ligand_name} – pose_idx is None")
            continue

        try:
            ligands = list(sdf_supplier(str(sdf_path)))
        except Exception as e:
            logger.warning(f"[ProLIF] Skipping {ligand_name} – failed to read SDF: {e}")
            continue

        if len(ligands) <= pose_idx:
            logger.warning(f"[ProLIF] Skipping {ligand_name} – pose_idx {pose_idx} out of range ({len(ligands)} poses)")
            continue

        mol = ligands[int(pose_idx)]
        if mol is None:
            logger.warning(f"[ProLIF] Skipping {ligand_name} – RDKit failed to parse SDF pose")
            continue

        if mol.GetNumConformers() == 0:
            logger.warning(f"[ProLIF] Skipping {ligand_name} – no conformers")
            continue

        # ---- ProLIF fingerprint ----
        try:
            ligand = plf.Molecule.from_rdkit(mol)
        except Exception as e:
            logger.warning(f"[ProLIF] Skipping {ligand_name} – failed to convert to ProLIF Molecule: {e}")
            continue

        try:
            fp = plf.Fingerprint()
            fp.run_from_iterable([ligand], protein, progress=False)
            df_fp_single = fp.to_dataframe()
        except Exception as e:
            logger.warning(f"[ProLIF] Skipping {ligand_name} – failed to compute fingerprint: {e}")
            continue

        # ---- Add to output ----
        dfs.append(df_fp_single)
        ligand_names.append(ligand_name)

    if not dfs:
        logger.warning("[ProLIF] No valid fingerprints generated; skipping ProLIF analysis.")
        return None

    # Concatenate into a single dataframe
    df_fp = pd.concat(dfs, axis=0)
    df_fp.index = ligand_names
    df_fp.index.name = "ligand"

    # Ensure boolean dtype, fill NaNs
    df_fp = df_fp.astype("boolean").fillna(False)

    return df_fp


def residue_number(res):
    """
    Extract residue number from labels like 'ARG250.' or 'THR52'
    """
    import re
    m = re.search(r"(\d+)", str(res))
    return int(m.group(1)) if m else float("inf")

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

    # Sort by residue number)
    residues = mat_plot.index.get_level_values(residue_level)

    unique_residues = sorted(
        set(residues),
        key=residue_number
    )

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
    ax.set_title("Barcode plot (prot x lig interactions)")

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
