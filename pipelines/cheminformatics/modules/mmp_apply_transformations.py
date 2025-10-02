from pipeline.task_registry import register_task
import json
import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.QED import qed
from pipeline.logger import setup_logger

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
    from rdkit.Chem import SaltRemover, rdmolops

    remover = SaltRemover.SaltRemover()
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

def is_valid_single_substituent_smarts(sub_smarts):
    try:
        mol = Chem.MolFromSmarts(sub_smarts)
        if mol is None:
            return False
        frags = Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False)
        return len(frags) == 1 and sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 0) == 1
    except:
        return False

MAX_VALENCE = 4  # Conservative cap for typical atoms

def get_dummy_neighbors(mol):
    """
    Find all dummy atoms and return a list of (dummy_idx, neighbor_idx) tuples.
    """
    dummy_pairs = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            neighbors = atom.GetNeighbors()
            if neighbors:
                dummy_pairs.append((atom.GetIdx(), neighbors[0].GetIdx()))
    return dummy_pairs


def reassemble_full_lead(lead_mol, core_smarts, sub_smarts):
    core_mol = Chem.MolFromSmarts(core_smarts)
    sub_mol = Chem.MolFromSmarts(sub_smarts)

    if core_mol is None or sub_mol is None:
        logger.debug(f"Invalid core or substituent SMARTS: {core_smarts}, {sub_smarts}")
        return None

    match = lead_mol.GetSubstructMatch(core_mol)
    if not match:
        logger.debug("No substructure match for core.")
        return None

    # === Find dummy atoms and their neighbors ===
    core_dummy_pairs = get_dummy_neighbors(core_mol)
    sub_dummy_pairs = get_dummy_neighbors(sub_mol)

    if len(core_dummy_pairs) != len(sub_dummy_pairs):
        logger.debug(f"Mismatch in dummy count: core={len(core_dummy_pairs)}, sub={len(sub_dummy_pairs)}")
        return None

    # Map core dummy indices to lead molecule indices
    lead_dummy_to_neighbor = []
    for dummy_idx_core, neighbor_idx_core in core_dummy_pairs:
        dummy_in_lead = match[dummy_idx_core]
        neighbor_in_lead = None
        for n in lead_mol.GetAtomWithIdx(dummy_in_lead).GetNeighbors():
            neighbor_in_lead = n.GetIdx()
            break
        if neighbor_in_lead is None:
            logger.debug("Dummy in lead has no neighbor.")
            return None
        lead_dummy_to_neighbor.append((dummy_in_lead, neighbor_in_lead))

    # === Remove dummy atoms from both lead and substituent ===
    lead_rw = Chem.RWMol(lead_mol)
    for dummy_idx, _ in sorted(lead_dummy_to_neighbor, reverse=True):
        lead_rw.RemoveAtom(dummy_idx)
    lead_clean = lead_rw.GetMol()

    sub_rw = Chem.RWMol(sub_mol)
    sub_dummy_to_neighbor = []
    for dummy_idx, neighbor_idx in sorted(sub_dummy_pairs, reverse=True):
        sub_dummy_to_neighbor.append((dummy_idx, neighbor_idx))
        sub_rw.RemoveAtom(dummy_idx)
    sub_clean = sub_rw.GetMol()

    # === Find updated neighbor indices in sub_clean ===
    updated_sub_neighbors = []
    old_to_new_idx = {}
    old_idx = 0
    for atom in sub_mol.GetAtoms():
        if atom.GetAtomicNum() != 0:
            old_to_new_idx[atom.GetIdx()] = old_idx
            old_idx += 1

    for _, old_neighbor_idx in sub_dummy_to_neighbor:
        new_idx = old_to_new_idx.get(old_neighbor_idx)
        if new_idx is None:
            logger.debug("Substituent neighbor mapping failed.")
            return None
        updated_sub_neighbors.append(new_idx)

    # === Valence checks ===
    for lead_neighbor_idx, sub_neighbor_idx in zip([n for _, n in lead_dummy_to_neighbor], updated_sub_neighbors):
        try:
            lead_atom = lead_clean.GetAtomWithIdx(lead_neighbor_idx)
            sub_atom = sub_clean.GetAtomWithIdx(sub_neighbor_idx)
            lead_valence = sum(b.GetBondTypeAsDouble() for b in lead_atom.GetBonds())
            sub_valence = sum(b.GetBondTypeAsDouble() for b in sub_atom.GetBonds())
            if lead_valence >= MAX_VALENCE or sub_valence >= MAX_VALENCE:
                logger.debug(f"Valence too high: lead={lead_valence}, sub={sub_valence}")
                return None
        except Exception as e:
            logger.debug(f"Valence check error: {e}")
            return None

    # === Combine and bond ===
    combo = Chem.CombineMols(lead_clean, sub_clean)
    em = Chem.EditableMol(combo)
    sub_offset = lead_clean.GetNumAtoms()

    for (lead_dummy_idx, lead_neighbor_idx), sub_neighbor_idx in zip(lead_dummy_to_neighbor, updated_sub_neighbors):
        sub_atom_idx = sub_offset + sub_neighbor_idx
        if lead_neighbor_idx >= combo.GetNumAtoms() or sub_atom_idx >= combo.GetNumAtoms():
            logger.debug(f"AddBond index out of range: {lead_neighbor_idx}, {sub_atom_idx}")
            return None
        try:
            em.AddBond(lead_neighbor_idx, sub_atom_idx, Chem.rdchem.BondType.SINGLE)
        except Exception as e:
            logger.debug(f"Failed to add bond: {e}")
            return None

    combined = em.GetMol()

    # === Sanitization ===
    try:
        Chem.SanitizeMol(combined, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
    except Exception as e:
        logger.debug(f"Partial sanitization failed: {e}")
        return None

    try:
        Chem.SanitizeMol(combined)
    except Exception as e:
        logger.debug(f"Full sanitization failed: {e}")
        return None

    # === Optional Filters ===
    if combined.GetNumAtoms() > 100:
        logger.debug("Filtered: molecule too large.")
        return None
    if '.' in Chem.MolToSmiles(combined):
        logger.debug("Filtered: disconnected molecule.")
        return None

    return combined


def score_core_mol(core_mol):
    num_atoms = core_mol.GetNumAtoms()
    num_bonds = core_mol.GetNumBonds()
    num_rings = core_mol.GetRingInfo().NumRings()
    num_wildcards = sum(1 for atom in core_mol.GetAtoms() if atom.GetAtomicNum() == 0)

    return (
        2.0 * num_atoms +
        1.5 * num_bonds +
        3.0 * num_rings -
        5.0 * num_wildcards
    )



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

    core_mol_map = {}
    for smarts, subs in scaffold_to_subs_raw.items():
        mol = Chem.MolFromSmarts(smarts)
        if mol is not None:
            try:
                Chem.GetSymmSSSR(mol)  # Initializes ring info
                core_mol_map[smarts] = {"mol": mol, "subs": subs}
            except Exception as e:
                logger.debug(f"Sanitization failed for core: {smarts} — {e}")
        else:
            logger.debug(f"Invalid SMARTS in rule file: {smarts}")

    df = pd.read_csv(transform_input_file)
    if "smiles" not in df.columns:
        raise ValueError(f"'smiles' column not found in {transform_input_file}")

    suggestions = []

    for smi in df["smiles"]:
        lead_mol = Chem.MolFromSmiles(smi)
        lead_mol = preprocess_molecule(lead_mol)
        if lead_mol is None:
            logger.debug(f"Skipping lead (preprocess failed): {smi}")
            continue

        # === Select the best-matching core ===
        best_core_smarts = None
        best_core_info = None
        best_score = -1

        for core_smarts, info in core_mol_map.items():
            core_mol = info["mol"]
            match = lead_mol.GetSubstructMatch(core_mol)
            if not match:
                continue

            # Prefer more complex cores: more bonds = more context
            score = score_core_mol(core_mol)

            if score > best_score:
                best_score = score
                best_core_smarts = core_smarts
                best_core_info = info

        if best_core_smarts is None:
            logger.debug(f"No matching core found for lead: {smi}")
            continue

        logger.debug(f"Best match found for core {best_core_smarts} in lead {smi}")

        for new_sub in best_core_info["subs"]:
            if not is_valid_single_substituent_smarts(new_sub):
                logger.debug(f"Skipping invalid substituent (not single molecule): {new_sub}")
                continue

            new_mol = reassemble_full_lead(lead_mol, best_core_smarts, new_sub)
            if new_mol is None:
                continue

            try:
                new_smi = Chem.MolToSmiles(new_mol, canonical=True)
                props = compute_physchem(new_mol)
                suggestions.append(props)
                logger.debug(f"{smi} => {new_smi} (core {best_core_smarts} → new_sub {new_sub})")
            except Exception as e:
                logger.debug(f"Failed to compute properties for transformed molecule: {e}")
                continue

    df_out = pd.DataFrame(suggestions)
    output_dir = Path(config.get("output", {}).get("directory", "outputs/mmp_based_transformation"))
    output_dir.mkdir(parents=True, exist_ok=True)

    if df_out.empty:
        logger.warning("No transformed molecules were generated.")
        return (Path(transform_input_file).stem, df_out)

    csv_path = output_dir / f"{Path(transform_input_file).stem}_mmp_based_transformation.csv"
    smi_path = output_dir / f"{Path(transform_input_file).stem}_mmp_based_transformation.smi"

    df_out.to_csv(csv_path, index=False)
    with open(smi_path, "w") as f:
        for smi_val in df_out["smiles"]:
            f.write(smi_val + "\n")

    logger.info(f"Generated {len(df_out)} transformed molecules from input: {transform_input_file}")
    logger.info(f"CSV written: {csv_path}")
    logger.info(f"SMILES file written: {smi_path}")

    return (Path(transform_input_file).stem, df_out)

