import subprocess
from pathlib import Path
from rdkit import Chem

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

import subprocess
from pathlib import Path

def run_fpocket(pdb_path: Path, output_dir: Path) -> Path:
    """
    Runs fpocket on the given PDB file and returns the path to the output directory.

    Args:
        pdb_path (Path): Path to the protein PDB file.
        output_dir (Path): Directory where fpocket output should live for downstream tasks.

    Returns:
        Path: Path to fpocket output directory (ends in *_out)
    """
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file does not exist: {pdb_path}")

    try:
        subprocess.run(["fpocket", "-f", str(pdb_path)], check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"fpocket failed to run on {pdb_path}") from e

    fpocket_dir = pdb_path.with_name(pdb_path.stem + "_out")

    if not fpocket_dir.exists():
        raise FileNotFoundError(f"Expected fpocket output dir not found: {fpocket_dir}")

    return fpocket_dir

def get_center_from_pdb(pdb_path: Path):
    xs, ys, zs = [], [], []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM", "pocket")):  # Include pocket if needed
                # PDB format fixed columns for x,y,z coords
                # Columns: x=30-37, y=38-45, z=46-53 (1-based indexing)
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    xs.append(x)
                    ys.append(y)
                    zs.append(z)
                except ValueError:
                    continue
    if not xs:
        raise RuntimeError(f"No atom coordinates found in {pdb_path}")
    center_x = sum(xs) / len(xs)
    center_y = sum(ys) / len(ys)
    center_z = sum(zs) / len(zs)
    return (center_x, center_y, center_z)

def extract_pocket_centers(fpocket_output_dir: Path, top_n: int):
    pockets_dir = fpocket_output_dir / "pockets"
    pocket_files = sorted(pockets_dir.glob("pocket*_atm.pdb"))

    # Limit to top N pockets (fpocket ranks pockets by score)
    selected_pockets = pocket_files[:top_n]

    centers = []
    for pocket_file in selected_pockets:
        # parse coords from pocket_file, for example average atom coords
        coords = get_center_from_pdb(pocket_file)  # You likely have this function
        centers.append(coords)

    return centers

def plot_multi_pocket_scores(csv_path: Path, output_dir: Path):
    """
    Plots violin plots of docking scores across pockets for each ligand.

    Args:
        csv_path (Path): Path to the summary CSV of docking scores.
        output_dir (Path): Where to save the plot.
    """
    df = pd.read_csv(csv_path)

    required_cols = {'ligand', 'score', 'pocket'}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"CSV must contain {required_cols} columns.")

    df["pocket"] = df["pocket"].astype(str)

    plt.figure(figsize=(12, 6))
    sns.violinplot(data=df, x="ligand", y="score", hue="pocket", inner=None, linewidth=1)
    sns.swarmplot(data=df, x="ligand", y="score", hue="pocket", dodge=True, size=3, color=".2")

    plt.title("Docking Affinities Across Pockets (minimizedAffinity)")
    plt.ylabel("Docking Score (kcal/mol)")
    plt.xticks(rotation=45)
    plt.tight_layout()

    handles, labels = plt.gca().get_legend_handles_labels()
    n = len(set(df["pocket"]))
    plt.legend(handles[:n], labels[:n], title="Pocket")

    plot_path = output_dir / "multi_pocket_scores.png"
    plt.savefig(plot_path)
    plt.close()