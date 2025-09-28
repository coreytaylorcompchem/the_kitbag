from pipeline.task_registry import register_task
from rdkit import Chem
from rdkit.Chem import rdMMPA, Descriptors, rdMolDescriptors
from rdkit.Chem.QED import qed
from pipeline.logger import setup_logger
from pathlib import Path
import pandas as pd
import json

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def compute_physchem(mol):
    return {
        "mol": mol,
        "smiles": Chem.MolToSmiles(mol, canonical=True),
        "mw": Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "hbd": Descriptors.NumHDonors(mol),
        "hba": Descriptors.NumHAcceptors(mol),
        "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
        "tpsa": Descriptors.TPSA(mol),
        "qed": qed(mol),
        "stereocenters": rdMolDescriptors.CalcNumAtomStereoCenters(mol),
    }

def rebuild_molecule(core_smiles, substituent_smiles):
    core_mol = Chem.MolFromSmiles(core_smiles)
    sub_mol = Chem.MolFromSmiles(substituent_smiles)
    if core_mol is None or sub_mol is None:
        return None

    # Find dummy atoms with isotope labels 1 in core, 2 in substituent
    core_attach = [a for a in core_mol.GetAtoms() if a.GetAtomicNum() == 0 and a.GetIsotope() == 1]
    sub_attach = [a for a in sub_mol.GetAtoms() if a.GetAtomicNum() == 0 and a.GetIsotope() == 2]
    if not core_attach or not sub_attach:
        return None

    combo = Chem.CombineMols(core_mol, sub_mol)
    emol = Chem.EditableMol(combo)

    core_idx = core_attach[0].GetIdx()
    sub_idx = sub_attach[0].GetIdx() + core_mol.GetNumAtoms()  # offset by core atoms count

    emol.AddBond(core_idx, sub_idx, Chem.rdchem.BondType.SINGLE)
    new_mol = emol.GetMol()

    # Remove dummy atoms in reverse order to keep indices valid
    dummy_indices = [atom.GetIdx() for atom in new_mol.GetAtoms() if atom.GetAtomicNum() == 0]
    rw_mol = Chem.RWMol(new_mol)
    for idx in sorted(dummy_indices, reverse=True):
        rw_mol.RemoveAtom(idx)
    new_mol = rw_mol.GetMol()

    try:
        Chem.SanitizeMol(new_mol)
    except Exception as e:
        logger.debug(f"Sanitization failed during rebuilding molecule: {e}")
        return None

    return new_mol

@register_task("mmp_apply_transformations", category="Generation", description="Applies learned MMP rules to lead molecules")
def mmp_apply_transformations(config, data=None):
    # Prefer input from mmp_apply_transformations config, fallback to top-level input_file
    transform_input_file = config.get("mmp_apply_transformations", {}).get("input_file") or config.get("input_file")
    if not transform_input_file or not Path(transform_input_file).exists():
        raise FileNotFoundError(f"Transformation input file not found: {transform_input_file}")

    rules_file = Path(config.get("output", {}).get("directory", "outputs/mmp_based_transformation")) / "mmp_rules.json"
    if not rules_file.exists():
        raise FileNotFoundError(f"Missing transformation rules file: {rules_file}")

    with open(rules_file, "r") as f:
        scaffold_to_subs = json.load(f)

    df = pd.read_csv(transform_input_file)
    if "smiles" not in df.columns:
        raise ValueError(f"'smiles' column not found in {transform_input_file}")

    suggestions = []

    for smi in df["smiles"]:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            logger.debug(f"Invalid molecule skipped: {smi}")
            continue

        try:
            # Fragment using the same settings as learning
            cuts = rdMMPA.FragmentMol(mol, maxCuts=2, maxCutBonds=60, pattern="[*:1]~[*:2]")
            if not cuts:
                # Optionally try fallback maxCuts=1
                cuts = rdMMPA.FragmentMol(mol, maxCuts=1, maxCutBonds=60, pattern="[*:1]~[*:2]")
                if not cuts:
                    logger.debug(f"No MMP cuts found for: {smi}")
                    continue

            for core_mol, sub_mol in cuts:
                if core_mol is None or sub_mol is None:
                    continue
                core_smiles = Chem.MolToSmiles(core_mol, canonical=True)
                sub_smiles = Chem.MolToSmiles(sub_mol, canonical=True)

                replacements = scaffold_to_subs.get(core_smiles, [])
                if not replacements:
                    logger.debug(f"Core {core_smiles} not found in learned rules.")
                    continue

                for new_sub in replacements:
                    if new_sub == sub_smiles:
                        continue
                    new_mol = rebuild_molecule(core_smiles, new_sub)
                    if new_mol is None:
                        continue
                    props = compute_physchem(new_mol)
                    suggestions.append(props)

        except Exception as e:
            logger.debug(f"Fragmentation or rebuilding failed for {smi}: {e}")
            continue

    df_out = pd.DataFrame(suggestions)
    logger.info(f"Generated {len(df_out)} transformed molecules from input: {transform_input_file}")

    return (Path(transform_input_file).stem, df_out)
