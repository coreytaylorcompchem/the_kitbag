from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdMolDescriptors, DataStructs
from rdkit.Chem import BRICS
from rdkit import RDLogger

import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

from pipeline.logger import setup_logger
from pipeline.task_registry import register_task

lg = RDLogger.logger()
lg.setLevel(RDLogger.ERROR)  # suppress RDKit warnings

logger = setup_logger(__name__, debug_mode=False, simple_format=True)


def get_murcko_scaffold(mol):
    try:
        return MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        logger.debug(f"Failed to get Murcko scaffold: {e}")
        return None


def mol_to_fingerprint(mol):
    # Returns ExplicitBitVect Morgan fingerprint (2048 bits, radius=2)
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)


def tanimoto(fp1, fp2):
    return DataStructs.TanimotoSimilarity(fp1, fp2)


def deduplicate_by_similarity(mols, threshold=0.7, max_to_keep=100):
    selected = []
    fps = []
    for mol in mols:
        fp = mol_to_fingerprint(mol)
        if all(tanimoto(fp, prev_fp) < threshold for prev_fp in fps):
            selected.append(mol)
            fps.append(fp)
        if len(selected) >= max_to_keep:
            break
    return selected


@register_task("fragment_novelty_filtering", category="Filtering", description="Scaffold-based fragment selection for novelty filtering")
def fragment_novelty_filtering(config, data=None):
    input_file = config.get("input_file")
    if input_file is None:
        raise ValueError("No input_file specified in config")

    input_path = Path(input_file)
    df = pd.read_csv(input_file)
    if "smiles" not in df.columns:
        raise ValueError("Input CSV must contain 'smiles' column.")

    scaffold_threshold = config.get("fragment_novelty_filtering", {}).get("scaffold_max_similarity", 0.7)
    max_fragments = config.get("fragment_novelty_filtering", {}).get("max_fragments_to_keep", 100)
    min_frag_atoms = config.get("fragment_novelty_filtering", {}).get("min_heavy_atoms", 8)

    scaffolds = []
    valid_mols = []

    for smi in df["smiles"].dropna():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        scaffold = get_murcko_scaffold(mol)
        if scaffold is None or scaffold.GetNumAtoms() < min_frag_atoms:
            continue
        scaffolds.append(scaffold)
        valid_mols.append(mol)

    if not scaffolds:
        logger.debug("No scaffolds meeting criteria were extracted from this chunk.")
        empty_df = pd.DataFrame(columns=["scaffold_smiles"])
        return (input_path.stem, empty_df)

    logger.debug(f"Extracted {len(scaffolds)} scaffolds before novelty filtering.")

    # Unique scaffold smiles before filtering for plotting
    unique_before = len(set(Chem.MolToSmiles(scaf, canonical=True) for scaf in scaffolds))

    novel_scaffolds = deduplicate_by_similarity(scaffolds, threshold=scaffold_threshold, max_to_keep=max_fragments)
    novel_smis = [Chem.MolToSmiles(scaf, canonical=True) for scaf in novel_scaffolds]
    unique_after = len(set(novel_smis))

    logger.debug(f"{unique_after} novel scaffolds retained after filtering.")

    # Plotting
    plt.figure(figsize=(6, 4))
    plt.bar(["Before Filtering", "After Filtering"], [unique_before, unique_after], color=["#4c72b0", "#55a868"])
    plt.ylabel("Unique Scaffold Count")
    plt.title("Scaffold Novelty Filtering")
    plt.tight_layout()

    output_dir = Path(config.get("output", {}).get("directory", "outputs"))
    plots_output_dir = output_dir / "_plots"
    
    output_dir.mkdir(parents=True, exist_ok=True)
    plots_output_dir.mkdir(parents=True, exist_ok=True) 

    plot_file = plots_output_dir / f"{input_path.stem}_scaffold_novelty.png"
    plt.savefig(plot_file)
    plt.close()
    logger.debug(f"Saved scaffold novelty plot: {plot_file}")

    # Save filtered scaffolds CSV as before
    result_df = pd.DataFrame({"scaffold_smiles": novel_smis})
    output_file = output_dir / config.get("output", {}).get("filename", "novel_fragments.csv")
    result_df.to_csv(output_file, index=False)

    logger.debug(f"Novel fragments written to: {output_file}")

    return (input_path.stem, result_df)
