from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdMolDescriptors, DataStructs
from rdkit.Chem import BRICS
from rdkit import RDLogger

lg = RDLogger.logger()
lg.setLevel(RDLogger.ERROR)  # or RDLogger.CRITICAL
import pandas as pd
from pathlib import Path

from pipeline.logger import setup_logger
from pipeline.task_registry import register_task

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

def fragment_molecule(mol):
    try:
        return list(BRICS.BRICSDecompose(mol, minFragmentSize=3))
    except Exception as e:
        logger.debug(f"BRICS decomposition failed: {e}")
        return []

@register_task("fragment_novelty_filtering", description="Scaffold-based fragment selection for novelty filtering")
def fragment_novelty_filtering(config, data=None):
    input_file = config.get("input_file")
    if input_file is None:
        raise ValueError("No input_file specified in config")

    input_path = Path(input_file)
    df = pd.read_csv(input_file)
    if "smiles" not in df.columns:
        raise ValueError("Input CSV must contain 'smiles' column.")

    scaffold_threshold = config.get("fragment_novelty_filtering", {}).get("scaffold_max_similarity", 0.7)
    max_fragments = config.get("fragment_novelty", {}).get("max_fragments_to_keep", 100)
    min_frag_atoms = config.get("fragment_novelty", {}).get("min_fragment_size", 8)

    valid_mols = []
    scaffolds = []

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
        logger.debug("No scaffolds meeting 'scaffold_max_similarity' threshold were extracted from this chunk.")
        return (input_path.stem, pd.DataFrame(columns=["smiles"]))

    logger.debug(f"Extracted {len(scaffolds)} scaffolds. Performing novelty filtering...")

    novel_scaffolds = deduplicate_by_similarity(scaffolds, threshold=scaffold_threshold, max_to_keep=max_fragments)
    novel_smis = [Chem.MolToSmiles(scaf, canonical=True) for scaf in novel_scaffolds]

    result_df = pd.DataFrame({"scaffold_smiles": novel_smis})

    output_dir = Path(config.get("output", {}).get("directory", "outputs"))
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / config.get("output", {}).get("filename", "novel_fragments.csv")

    result_df.to_csv(output_file, index=False)

    logger.debug(f"Novel fragments written to: {output_file}")

    return (input_path.stem, result_df)

