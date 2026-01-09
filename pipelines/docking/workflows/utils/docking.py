import csv
import re

from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
from collections import OrderedDict

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Rectangle

from rdkit import Chem
from rdkit.Chem.rdmolfiles import MolFromPDBFile

import MDAnalysis as mda

import prolif as plf
from prolif import sdf_supplier

from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

def generate_ligands_csv_from_txt(txt_path: Path, csv_path: Path):
    with open(txt_path, 'r') as f:
        smiles_list = [line.strip() for line in f if line.strip()]

    if not smiles_list:
        raise ValueError(f"No SMILES found in {txt_path}")

    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['name', 'smiles'])
        for i, smi in enumerate(smiles_list, start=1):
            writer.writerow([f"ligand{i}", smi])

    logger.info(f"Generated ligands.csv at {csv_path} with {len(smiles_list)} ligands.")

def validate_ligands_csv(csv_path: Path):
    if not csv_path.exists():
        raise FileNotFoundError(f"Expected ligands.csv at {csv_path} not found.")
    with open(csv_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        if 'name' not in reader.fieldnames or 'smiles' not in reader.fieldnames:
            raise ValueError(f"ligands.csv must have 'name' and 'smiles' columns.")

def get_docking_box(config, protein_preparer):
    """
    Determine the docking box center and size from config.
    """
    docking_cfg = config.get('docking', {})

    if 'center' in docking_cfg:
        center = tuple(docking_cfg['center'])

    elif 'ref_ligand_path' in docking_cfg:
        ref_path = docking_cfg['ref_ligand_path']
        ligand_mol = Chem.MolFromMolFile(ref_path)
        if ligand_mol is None:
            raise ValueError(f"[ERROR] Reference ligand file invalid: {ref_path}")
        conf = ligand_mol.GetConformer()
        pts = [list(conf.GetAtomPosition(i)) for i in range(ligand_mol.GetNumAtoms())]
        center = tuple(sum(x) / len(x) for x in zip(*pts))

    elif docking_cfg.get('use_crystal_ligand', False):
        ref_lig_path = protein_preparer.reference_ligand_path
        ligand_mol = Chem.MolFromPDBFile(str(ref_lig_path))
        if ligand_mol is None:
            raise ValueError("[ERROR] Failed to read crystallized ligand from protein.")
        conf = ligand_mol.GetConformer()
        pts = [list(conf.GetAtomPosition(i)) for i in range(ligand_mol.GetNumAtoms())]
        center = tuple(sum(x) / len(x) for x in zip(*pts))

    else:
        raise ValueError(
            "[ERROR] Docking center not specified. Provide one of:\n"
            "  - center: [x, y, z]\n"
            "  - ref_ligand_path: /path/to/ligand\n"
            "  - use_crystal_ligand: true"
        )

    if "size" not in docking_cfg:
        raise ValueError("[ERROR] Docking 'size' must be specified in config['docking'].")

    size = tuple(docking_cfg["size"])

    return center, size

def validate_config(config, required_fields):
    missing_fields = []
    for field in required_fields:
        parts = field.split('.')
        value = config
        for part in parts:
            if isinstance(value, dict) and part in value:
                value = value[part]
            else:
                value = None
                break
        if value in (None, [], ''):
            missing_fields.append(field)
    
    if missing_fields:
        logger.warning(f"Required config fields are missing: {missing_fields}")

def summarise_docking_results(output_dir: Path, summary_csv_path: Path):
    """
    Aggregates docking scores from per-pocket result files into one summary CSV.

    Works for both flat and nested pocket structures.

    Args:
        output_dir (Path): Root output directory.
        summary_csv_path (Path): Where to save the combined CSV.
    """
    all_rows = []

    # Recursively search for all docking_scores.csv files
    for pocket_dir in output_dir.rglob("pocket_*"):
        if not pocket_dir.is_dir():
            continue

        docking_results_file = pocket_dir / "docking_scores.csv"
        if docking_results_file.exists():
            try:
                df = pd.read_csv(docking_results_file)

                pocket_id = pocket_dir.name.replace("pocket_", "")
                df['pocket'] = pocket_id

                # If nested, grab parent folder name as structure ID
                parent = pocket_dir.parent
                if parent != output_dir:
                    df['structure'] = parent.name
                else:
                    df['structure'] = None

                all_rows.append(df)

            except Exception as e:
                logger.warning(f"Failed to read {docking_results_file}: {e}")

    if all_rows:
        combined = pd.concat(all_rows, ignore_index=True)
        combined.to_csv(summary_csv_path, index=False)
        logger.info(f"Docking results summarised into: {summary_csv_path}")
    else:
        logger.info("No docking result files found.")

def extract_crystal_ligand_center(pdb_path, ligand_resname=None):
    """
    Extracts the 3D center of the co-crystallized ligand.
    Assumes it's the only non-protein residue or that ligand_resname is provided.
    """
    mol = MolFromPDBFile(str(pdb_path), removeHs=False)
    if mol is None:
        raise ValueError(f"Could not read molecule from {pdb_path}")

    conformer = mol.GetConformer()
    atom_positions = []
    for atom in mol.GetAtoms():
        resinfo = atom.GetPDBResidueInfo()
        if resinfo is None:
            continue
        if ligand_resname and resinfo.GetResidueName().strip() != ligand_resname:
            continue
        if resinfo.GetResidueName().strip() in ["HOH", "WAT"]:
            continue
        atom_positions.append(conformer.GetAtomPosition(atom.GetIdx()))

    if not atom_positions:
        raise ValueError("No ligand atoms found to determine center.")

    coords = np.array([[pos.x, pos.y, pos.z] for pos in atom_positions])
    center = coords.mean(axis=0)
    return center.tolist()

def plot_multi_structure_scores(summary_csv_path: Path, output_dir: Path):
    """
    Plots ligand performance across multiple structures/pockets.

    - Heatmap of best scores
    - Boxplot of score distributions
    """
    df = pd.read_csv(summary_csv_path)

    if 'structure' not in df.columns:
        logger.warning("Missing 'structure' column. Skipping plots.")
        return

    output_dir.mkdir(exist_ok=True, parents=True)

    # drop missing scores
    df = df.dropna(subset=["score"])

    heatmap_df = df.groupby(['ligand', 'structure'])['score'].min().unstack()

    plt.figure(figsize=(12, 8))
    sns.heatmap(heatmap_df, annot=True, cmap="coolwarm_r", fmt=".2f")
    plt.title("Best Docking Score per Ligand vs Structure")
    plt.ylabel("Ligand")
    plt.xlabel("Structure")
    plt.tight_layout()
    plt.savefig(output_dir / "ligand_structure_heatmap.png", dpi=300)
    plt.close()

    plt.figure(figsize=(14, 6))
    sns.boxplot(data=df, x='ligand', y='score', hue='structure')
    plt.title("Docking Score Distribution per Ligand")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_dir / "ligand_score_boxplot.png", dpi=300)
    plt.close()

    logger.info("Generated ligand-structure docking plots.")

def extract_scores_from_docking_output(output_path):
    """
    Dispatch score extraction based on file extension.
    Accepts a single Path or a list of Paths.
    Returns a list of score dictionaries.
    """
    # Unwrap list if needed
    if isinstance(output_path, list):
        if not output_path:
            logger.warning("Empty docking output path list.")
            return []
        output_path = output_path[0]

    # Ensure Path object
    if isinstance(output_path, str):
        output_path = Path(output_path)

    if not isinstance(output_path, Path):
        logger.warning(f"Invalid docking output type: {type(output_path)}")
        return []

    ext = output_path.suffix.lower()
    if ext == ".sdf":
        return extract_scores_from_gnina_sdf(output_path)
    elif ext == ".pdbqt":
        return extract_scores_from_unidock_pdbqt(output_path)
    else:
        logger.warning(f"Unknown docking output format: {output_path}")
        return []

def extract_scores_from_gnina_sdf(sdf_path: Path):
    """
    Extract docking scores from a GNINA SDF output file.
    """
    results = []
    try:
        suppl = Chem.SDMolSupplier(str(sdf_path))
    except Exception as e:
        logger.error(f"Failed to read SDF file {sdf_path}: {e}")
        return results

    if not suppl:
        logger.warning(f"RDKit could not load molecules from SDF: {sdf_path}")
        return results

    for pose_idx, mol in enumerate(suppl):
        if mol is None:
            logger.warning(f"Skipping invalid molecule at pose {pose_idx} in {sdf_path}")
            continue
        try:
            score = float(mol.GetProp("minimizedAffinity"))
            cnn_score = float(mol.GetProp("CNNscore"))
        except KeyError:
            logger.warning(f"Missing scores in pose {pose_idx} of {sdf_path}")
            continue
        except Exception as e:
            logger.warning(f"Error parsing scores in pose {pose_idx} of {sdf_path}: {e}")
            continue

        results.append({
            "score": score,
            "pose_rank": cnn_score,
            "pose_idx": pose_idx
        })
    return results

def extract_scores_from_unidock_pdbqt(pdbqt_path: Path):
    """
    Extract scores from Uni-Dock PDBQT output by parsing REMARK lines.
    Assumes docking poses are separated by `MODEL` / `ENDMDL`.
    """
    logger.debug(f"Parsing Uni-Dock output: {pdbqt_path}")
    results = []
    pose_idx = -1

    with open(pdbqt_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("MODEL"):
                pose_idx += 1

            elif line.startswith("REMARK") and "VINA" in line.upper() and "RESULT" in line.upper():
                try:
                    # Extract all float-like numbers from the line
                    parts = line.replace(":", " ").split()
                    float_vals = [float(p) for p in parts if _is_float(p)]

                    if float_vals:
                        score = float_vals[0]  # The binding affinity is usually the first
                        results.append({
                            "score": score,
                            "pose_rank": None,  # Uni-Dock may not output CNNscore
                            "pose_idx": pose_idx
                        })
                except Exception as e:
                    logger.warning(f"Failed to parse score in {pdbqt_path} at pose {pose_idx}: {e}")
    return results

def _is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
    
# plotting code

def ligand_sort_key(name):
    """
    Extract integer index from ligand names like 'ligand12'.
    Falls back to name if no integer found.
    """
    m = re.search(r"(\d+)$", str(name))
    return int(m.group(1)) if m else name


def plot_docking_scatter_with_errorbars(csv_path: Path, output_dir: Path):
    df = pd.read_csv(csv_path)

    # enforce ligand order by name (ligand1, ligand2, ...)
    df["ligand_order"] = df["name"].apply(ligand_sort_key)
    df = df.sort_values("ligand_order")

    # keep only valid docking scores
    df = df.dropna(subset=["docking_score"])
    df = df[df["docking_score"] < 0]  # remove invalid positive scores

    # aggregate per ligand
    agg = (
        df.groupby("name")
        .agg(
            min_score=("docking_score", "min"),
            max_score=("docking_score", "max"),
            best_score=("docking_score", "min"),  # best = most negative
        )
        .reset_index()
    )

    # assign ligand index (NOT docking-ranked)
    agg["mol_idx"] = agg["name"].apply(ligand_sort_key)

    # sort by ligand index
    agg = agg.sort_values("mol_idx")

    # correct asymmetric error bars
    lower_err = agg["best_score"] - agg["min_score"]
    upper_err = agg["max_score"] - agg["best_score"]
    yerr = [lower_err.abs(), upper_err.abs()]

    plt.figure(figsize=(10, 6))
    plt.errorbar(
        agg["mol_idx"],
        agg["best_score"],
        yerr=yerr,
        fmt="o",
        ecolor="gray",
        elinewidth=1,
        capsize=3,
        markersize=5
    )

    # reference lines
    plt.axhline(-7.0, linestyle="--", color="red", alpha=0.7, label="Weak binder")
    plt.axhline(-8.0, linestyle="--", color="orange", alpha=0.7, label="Good binder")
    plt.axhline(-9.0, linestyle="--", color="green", alpha=0.7, label="Strong binder")

    plt.xlabel("Molecule (ligand index)")
    plt.ylabel("Docking score (kcal/mol)")
    plt.title("Docking score range per ligand")

    ax = plt.gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    plt.legend()
    plt.tight_layout()

    out = output_dir / "docking_scatter_errorbars.png"
    plt.savefig(out, dpi=300)
    plt.close()


ADME_THRESHOLDS = {
    "cyp3a4": 7.0,
    "caco2": 20.0,
    "logd": 3.0,
    "herg": 0.5
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

            ax.axhline(
                threshold,
                linestyle="--",
                color="red",
                linewidth=1
            )

            ax.set_title(f"{adme.upper()} coloured by {desc}")
            ax.set_ylabel(adme)

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

# Interaction analysis

def select_best_poses(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)

    df = df.dropna(subset=["docking_score"])
    df = df[df["docking_score"] < 0]

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

    for _, row in best_pose_df.iterrows():
        sdf_path = row.get("sdf_path")
        pose_idx = row.get("pose_idx")

        if not sdf_path or not Path(sdf_path).exists():
            continue
        if pose_idx is None:
            continue

        ligands = list(sdf_supplier(str(sdf_path)))
        if pose_idx >= len(ligands):
            continue

        ligand = plf.Molecule.from_rdkit(ligands[int(pose_idx)])

        fp = plf.Fingerprint()
        fp.run_from_iterable([ligand], protein)

        df_fp_single = fp.to_dataframe()
        dfs.append(df_fp_single)
        ligand_names.append(row["name"])

    if not dfs:
        raise ValueError("No valid ProLIF fingerprints generated.")

    # Build DataFrame directly avoiding Prolif's crap helpers

    df_fp = pd.concat(dfs, axis=0)
    df_fp.index = ligand_names
    df_fp.index.name = "ligand"
    df_fp = df_fp.fillna(False)

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
        mat.groupby(level=interaction_level, axis=1)
           .mean()
           .mean(axis=0)
    )

    keep_interactions = interaction_freq[
        interaction_freq >= min_frequency
    ].index

    mat = mat.loc[
        :,
        mat.columns.get_level_values(interaction_level)
                   .isin(keep_interactions)
    ]

    mat = mat.sort_index(axis=1)
    mat_plot = mat.T  # residues × ligands


    # Build residue groups

    residues = mat_plot.index.get_level_values(residue_level)
    unique_residues = list(OrderedDict.fromkeys(residues))

    y_positions = []
    y_labels = []
    y_ticks = []

    y = 0
    row_map = []

    for res in unique_residues:
        idxs = [
            i for i, r in enumerate(residues) if r == res
        ]
        start = y
        for i in idxs:
            row_map.append((y, i))
            y += 1
        end = y - 1

        y_ticks.append((start + end) / 2)
        y_labels.append(res)

    # Plotting with matplotlib

    fig_w = max(8, mat_plot.shape[1] * 0.6)
    fig_h = max(6, len(row_map) * 0.25)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # background
    ax.imshow(
        np.zeros((len(row_map), mat_plot.shape[1])),
        cmap="Greys",
        vmin=0,
        vmax=1,
        aspect="auto"
    )

    # draw interactions
    for y_pos, row_idx in row_map:
        _, residue, interaction = mat_plot.index[row_idx]
        color = interaction_colors.get(interaction, "#000000")

        for x, ligand in enumerate(mat_plot.columns):
            if mat_plot.iloc[row_idx, x]:
                ax.add_patch(
                    Rectangle(
                        (x - 0.5, y_pos - 0.5),
                        1,
                        1,
                        facecolor=color,
                        edgecolor="none"
                    )
                )

    # FOrmatting axes similar to Prolif

    ax.set_xticks(np.arange(mat_plot.shape[1]))
    ax.set_xticklabels(mat_plot.columns, rotation=45, ha="right")

    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels)

    ax.set_xlabel("Ligand ID")
    ax.set_ylabel("Protein residue")

    ax.set_title("Barcode (P-L interactions)")

    # grid
    ax.set_xticks(np.arange(-.5, mat_plot.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(row_map), 1), minor=True)
    ax.grid(which="minor", color="lightgrey", linestyle="-", linewidth=0.5)


    # Build legend

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