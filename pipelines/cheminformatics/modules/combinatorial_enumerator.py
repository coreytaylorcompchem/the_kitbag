from rdkit import Chem
from rdkit.Chem import BRICS
import pandas as pd
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger
# from workflows import register_workflow

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_task("combinatorial_enumerator", category="Library generation", description="Generate combinatorial library from fragments and R-groups (WIP).")
def combinatorial_enumerator(config: dict, data: dict = None) -> dict:
    if data is None or "df" not in data:
        raise ValueError("Input data must contain a DataFrame under 'df'")
    
    df = data["df"]
    if "smiles" not in df.columns:
        raise ValueError("Input DataFrame must contain a 'smiles' column")

    max_molecules = config.get("max_molecules", 1000)
    min_fragments = config.get("min_fragments", 2)
    top_fragments = config.get("top_fragments", 20)

    logger.info(f"[combinatorial_enumerator] Extracting fragments from {len(df)} molecules...")

    mols = [Chem.MolFromSmiles(smi) for smi in df["smiles"].dropna()]
    brics_fragments = set()

    for mol in mols:
        try:
            frags = BRICS.BRICSDecompose(mol)
            brics_fragments.update(frags)
        except Exception as e:
            logger.debug(f"BRICS decomposition failed: {e}")

    logger.info(f"[combinatorial_enumerator] Found {len(brics_fragments)} unique BRICS fragments")

    if len(brics_fragments) < min_fragments:
        logger.warning("[combinatorial_enumerator] Not enough fragments for enumeration")
        return {"df": pd.DataFrame(), "mols": []}

    logger.info(f"[combinatorial_enumerator] Enumerating up to {max_molecules} new molecules")
    new_mols = []
    new_smiles = []

    try:
        combos = BRICS.BRICSBuild(list(brics_fragments))
        for i, mol in enumerate(combos):
            if i >= max_molecules:
                break
            if mol is None:
                continue
            smi = Chem.MolToSmiles(mol)
            new_mols.append(mol)
            new_smiles.append(smi)
    except Exception as e:
        logger.error(f"BRICS enumeration failed: {e}")
        return {"df": pd.DataFrame(), "mols": []}

    logger.info(f"[combinatorial_enumerator] Generated {len(new_smiles)} new molecules")

    new_df = pd.DataFrame({"smiles": new_smiles})
    return {"df": new_df, "mols": new_mols}
