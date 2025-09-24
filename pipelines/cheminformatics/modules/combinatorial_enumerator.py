from rdkit import Chem
from rdkit.Chem import BRICS
import pandas as pd
import gc
from collections import Counter
from tqdm import tqdm 
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_task("combinatorial_enumerator", category="Library generation", description="Generate combinatorial library from fragments and R-groups.")
def combinatorial_enumerator(config: dict, data: dict = None) -> dict:

    # Read combinatorial-specific parameters from yaml
    params = config.get("combinatorial_enumerator", {})

    if data is None or "df" not in data:
        raise ValueError("Input data must contain a DataFrame under 'df'")

    df = data["df"]
    if "smiles" not in df.columns:
        raise ValueError("Input DataFrame must contain a 'smiles' column")

    max_molecules = params.get("max_molecules", 100)
    min_fragments = params.get("min_fragments", 2)
    top_fragments = params.get("top_fragments", 20)
    max_attachment_points = params.get("max_attachment_points", 4)
    top_scaffolds = params.get("top_scaffolds", 100)
    top_rgroups = params.get("top_rgroups", 100)

    logger.debug(
        f"Parameters: max_molecules={max_molecules}, "
        f"min_fragments={min_fragments}, top_fragments={top_fragments}, "
        f"max_attachment_points={max_attachment_points}, top_scaffolds={top_scaffolds}, top_rgroups={top_rgroups}"
    )
    logger.info(f"Extracting fragments from {len(df)} molecules...")

    mols = [Chem.MolFromSmiles(smi) for smi in df["smiles"].dropna()]
    fragment_counter = Counter()

    for mol in tqdm(mols, desc="Fragmenting molecules", unit="mol"):
        try:
            frags = BRICS.BRICSDecompose(mol)
            fragment_counter.update(frags)
        except Exception as e:
            logger.debug(f"BRICS decomposition failed: {e}")

    if not fragment_counter:
        logger.warning("No fragments extracted.")
        return {"df": pd.DataFrame(), "mols": []}

    # Filter top fragments
    top_fragments_list = [frag for frag, _ in fragment_counter.most_common(top_fragments)]

    # Split fragments into scaffolds and R-groups
    scaffolds = [f for f in top_fragments_list if f.count('*') >= 2]
    r_groups = [f for f in top_fragments_list if f.count('*') == 1]

    # Filter based on max attachment points
    scaffolds = [f for f in scaffolds if f.count('*') <= max_attachment_points]

    # Limit number of scaffolds and R-groups
    scaffolds = scaffolds[:top_scaffolds]
    r_groups = r_groups[:top_rgroups]

    logger.info(f"Selected {len(scaffolds)} scaffolds and {len(r_groups)} R-groups")

    if len(scaffolds) < 1 or len(r_groups) < 1:
        logger.warning("Not enough fragments for enumeration")
        return {"df": pd.DataFrame(), "mols": []}

    logger.info(f"Enumerating up to {max_molecules} new molecules...")

    new_mols = []
    new_smiles = []

    # Combine all fragments for BRICSBuild
    build_frags = scaffolds + r_groups

    def limited_combos(fragments, max_mols):
        for i, mol in enumerate(BRICS.BRICSBuild(fragments)):
            if i >= max_mols:
                break
            if mol is not None:
                yield mol

    try:
        fragment_mols = [
            mol for mol in (Chem.MolFromSmiles(frag) for frag in build_frags)
            if mol is not None
        ]

        # Free memory
        del fragment_counter
        del top_fragments_list
        del mols
        gc.collect()

        for mol in tqdm(limited_combos(fragment_mols, max_molecules), desc="Enumerating molecules", unit="mol"):
            smi = Chem.MolToSmiles(mol)
            new_mols.append(mol)
            new_smiles.append(smi)

        del fragment_mols
        gc.collect()
    except Exception as e:
        logger.error(f"BRICS enumeration failed: {e}")
        return {"df": pd.DataFrame(), "mols": []}

    logger.info(f"Generated {len(new_smiles)} new molecules")

    new_df = pd.DataFrame({"smiles": new_smiles})
        
    return {"df": new_df, "mols": new_mols}
