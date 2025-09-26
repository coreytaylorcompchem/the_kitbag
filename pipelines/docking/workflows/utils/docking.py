import csv
from pathlib import Path
from rdkit import Chem

from rdkit.Chem.rdmolfiles import MolFromPDBFile
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

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