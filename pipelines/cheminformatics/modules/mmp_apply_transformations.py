from pipeline.task_registry import register_task
from rdkit import Chem
from rdkit.Chem import BRICS, Descriptors, rdMolDescriptors
from rdkit.Chem.QED import qed
from pipeline.logger import setup_logger
from pathlib import Path
import pandas as pd
import json

logger = setup_logger(__name__, debug_mode=True, simple_format=True)

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

def preprocess_molecule(mol):
    from rdkit.Chem import SaltRemover
    remover = SaltRemover.SaltRemover()
    from rdkit.Chem import rdmolops
    if mol is None:
        return None
    frags = rdmolops.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return None
    mol = max(frags, key=lambda m: m.GetNumAtoms())
    mol = remover.StripMol(mol, dontRemoveEverything=True)
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        try:
            mol.UpdatePropertyCache(strict=False)
            Chem.GetSymmSSSR(mol)
        except:
            return None
    return mol

def normalize_dummy_atoms(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return smi
    rw = Chem.RWMol(mol)
    dummy_atoms = sorted([a for a in rw.GetAtoms() if a.GetAtomicNum() == 0], key=lambda a: a.GetIdx())
    for i, atom in enumerate(dummy_atoms, start=1):
        atom.SetIsotope(i)
    return Chem.MolToSmiles(rw, canonical=True)

def reassemble_from_brics(core_smi, sub_smi):
    core = Chem.MolFromSmiles(core_smi)
    sub = Chem.MolFromSmiles(sub_smi)
    if core is None or sub is None:
        logger.debug(f"Cannot parse core or sub: {core_smi}, {sub_smi}")
        return None

    core_dummies = [a for a in core.GetAtoms() if a.GetAtomicNum() == 0]
    sub_dummies = [a for a in sub.GetAtoms() if a.GetAtomicNum() == 0]
    if not core_dummies or not sub_dummies:
        logger.debug(f"No dummy atoms in core or sub: {core_smi}, {sub_smi}")
        return None

    combo = Chem.CombineMols(core, sub)
    edcombo = Chem.EditableMol(combo)

    core_idx = core_dummies[0].GetIdx()
    sub_idx = sub_dummies[0].GetIdx() + core.GetNumAtoms()
    edcombo.AddBond(core_idx, sub_idx, Chem.BondType.SINGLE)

    newmol = edcombo.GetMol()

    dummy_idxs = [a.GetIdx() for a in newmol.GetAtoms() if a.GetAtomicNum() == 0]
    rw = Chem.RWMol(newmol)
    for idx in sorted(dummy_idxs, reverse=True):
        rw.RemoveAtom(idx)
    newmol = rw.GetMol()

    try:
        Chem.SanitizeMol(newmol)
    except Exception as e:
        logger.debug(f"Sanitize failed for {core_smi} + {sub_smi}: {e}")
        return None

    return newmol

def split_fragment_to_single_dummy(mol):
    """
    Given a fragment mol that may have multiple dummy atoms, 
    split it into subfragments each containing exactly one dummy atom.
    Returns list of SMILES strings.
    """
    dummy_idxs = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
    if len(dummy_idxs) <= 1:
        return [Chem.MolToSmiles(mol, canonical=True)]
    
    # Remove bonds between dummy atoms to split the fragment
    rw = Chem.RWMol(mol)
    bonds_to_break = []
    for bond in mol.GetBonds():
        begin = bond.GetBeginAtom().GetIdx()
        end = bond.GetEndAtom().GetIdx()
        if (begin in dummy_idxs) and (end in dummy_idxs):
            bonds_to_break.append(bond.GetIdx())
    # Remove bonds between dummy atoms
    for bidx in sorted(bonds_to_break, reverse=True):
        rw.RemoveBond(rw.GetBondWithIdx(bidx).GetBeginAtomIdx(), rw.GetBondWithIdx(bidx).GetEndAtomIdx())
    split_mol = rw.GetMol()

    # Get connected components
    frags = Chem.GetMolFrags(split_mol, asMols=True, sanitizeFrags=True)
    # Filter fragments with exactly one dummy atom
    single_dummy_frags = []
    for f in frags:
        dummies = [a for a in f.GetAtoms() if a.GetAtomicNum() == 0]
        if len(dummies) == 1:
            single_dummy_frags.append(Chem.MolToSmiles(f, canonical=True))
    return single_dummy_frags

@register_task("mmp_apply_transformations", category="Generation", description="Apply BRICS‑derived transformations to leads")
def mmp_apply_transformations(config, data=None):
    transform_input_file = config.get("mmp_apply_transformations", {}).get("input_file")
    if not transform_input_file or not Path(transform_input_file).exists():
        raise FileNotFoundError(f"Transformation input file not found: {transform_input_file}")

    rules_file = Path(config.get("output", {}).get("directory", "outputs/mmp_based_transformation")) / "mmp_rules_brics.json"
    if not rules_file.exists():
        raise FileNotFoundError(f"Missing BRICS transformation rules file: {rules_file}")

    with open(rules_file, "r") as f:
        scaffold_to_subs_raw = json.load(f)

    # Helper to get canonical smiles or fallback to original
    def normalize_smi(smi):
        mol = Chem.MolFromSmiles(smi)
        if mol:
            return Chem.MolToSmiles(mol, canonical=True)
        else:
            logger.debug(f"Failed to parse smiles for normalization: {smi}")
            return smi

    def count_heavy_atoms(smi):
        mol = Chem.MolFromSmiles(smi)
        return mol.GetNumHeavyAtoms() if mol else 0

    # Normalize rules keys and subs; filter out small cores (<6 heavy atoms)
    scaffold_to_subs = {}
    for core_smi_raw, subs_raw in scaffold_to_subs_raw.items():
        core_smi = normalize_smi(core_smi_raw)
        # if count_heavy_atoms(core_smi) < 6:
        #     logger.debug(f"Skipping small core in rules: {core_smi}")
        #     continue
        normalized_subs = set()
        for sub_smi_raw in subs_raw:
            sub_smi = normalize_smi(sub_smi_raw)
            normalized_subs.add(sub_smi)
        scaffold_to_subs[core_smi] = normalized_subs

    # Prebuild core mols map from normalized rules
    core_mol_map = {}
    for core_smi in scaffold_to_subs.keys():
        mol = Chem.MolFromSmiles(core_smi)
        if mol:
            core_mol_map[core_smi] = mol
        else:
            logger.debug(f"Invalid core smiles in learned rules: {core_smi}")

    df = pd.read_csv(transform_input_file)
    if "smiles" not in df.columns:
        raise ValueError(f"'smiles' column not found in {transform_input_file}")

    suggestions = []

    for smi in df["smiles"]:
        mol = Chem.MolFromSmiles(smi)
        mol = preprocess_molecule(mol)
        if mol is None:
            logger.debug(f"Skipping invalid lead: {smi}")
            continue

        try:
            brics_frags_raw = list(BRICS.BRICSDecompose(mol))
            # Normalize and filter lead fragments for cores (>=6 heavy atoms)
            brics_frags = [normalize_dummy_atoms(frag) for frag in brics_frags_raw]
            logger.debug(f"BRICS fragments for lead {smi} (normalized): {brics_frags}")
            # brics_frags = [f for f in brics_frags if count_heavy_atoms(f) >= 6]

            # logger.debug(f"BRICS fragments for lead {smi} (filtered): {brics_frags}")

            if not brics_frags:
                logger.debug(f"No BRICS fragments for lead: {smi}")
                continue

            # Map lead fragments to learned cores using strict isomorphism check (both directions)
            frag_to_core_map = {}
            possible_cores = set()

            for frag_smi in brics_frags:
                frag_mol = Chem.MolFromSmiles(frag_smi)
                if frag_mol is None:
                    logger.debug(f"Skipping invalid fragment in lead: {frag_smi}")
                    continue

                matched_core = None
                for core_smi, core_mol in core_mol_map.items():
                    # Match both ways to ensure exact core-fragment equivalence
                    if frag_mol.HasSubstructMatch(core_mol) and core_mol.HasSubstructMatch(frag_mol):
                        matched_core = core_smi
                        break

                if matched_core:
                    frag_to_core_map[frag_smi] = matched_core
                    possible_cores.add(matched_core)
                else:
                    logger.debug(f"Fragment {frag_smi} did not match any learned core")

            if not possible_cores:
                logger.debug(f"No matching cores in learned rules for lead: {smi}")
                continue

            # For each core matched in the lead, attempt substituent replacements
            for core_smi in possible_cores:
                learned_subs = scaffold_to_subs.get(core_smi, set())
                if not learned_subs:
                    continue

                # Consider all fragments in lead that match this core to find substituents
                for orig_sub in brics_frags:
                    if orig_sub == core_smi:
                        continue
                    # Only process substituents whose matched core matches current core_smi
                    if frag_to_core_map.get(orig_sub) != core_smi:
                        continue
                    if orig_sub not in learned_subs:
                        logger.debug(f"Original substituent {orig_sub} not in learned substitutions for core {core_smi}")
                        continue

                    for new_sub in learned_subs:
                        if new_sub == orig_sub:
                            continue

                        newmol = reassemble_from_brics(core_smi, normalize_dummy_atoms(new_sub))
                        if newmol is None:
                            continue

                        smi_new = Chem.MolToSmiles(newmol, canonical=True)
                        props = compute_physchem(newmol)
                        suggestions.append(props)

                        logger.debug(f"{smi}  => {smi_new} (core {core_smi}, replaced {orig_sub} → {new_sub})")

        except Exception as e:
            logger.debug(f"Error processing lead {smi}: {e}")
            continue

    df_out = pd.DataFrame(suggestions)

    output_dir = Path(config.get("output", {}).get("directory", "outputs/mmp_based_transformation"))
    output_dir.mkdir(parents=True, exist_ok=True)

    if df_out.empty:
        logger.warning("No transformed molecules were generated.")
        return (Path(transform_input_file).stem, df_out)

    csv_path = output_dir / f"{Path(transform_input_file).stem}_mmp_based_transformation.csv"
    df_out.to_csv(csv_path, index=False)

    smi_path = output_dir / f"{Path(transform_input_file).stem}_mmp_based_transformation.smi"
    with open(smi_path, "w") as f:
        for smi_val in df_out["smiles"]:
            f.write(smi_val + "\n")

    logger.info(f"Generated {len(df_out)} transformed molecules from input: {transform_input_file}")
    logger.info(f"CSV written: {csv_path}")
    logger.info(f"SMILES file written: {smi_path}")

    return (Path(transform_input_file).stem, df_out)



