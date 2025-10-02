from pipeline.task_registry import register_task
import math, json
from pathlib import Path
import pandas as pd
from tqdm import tqdm
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor

from rdkit import Chem
from rdkit.Chem import BRICS, rdmolops, SaltRemover
from rdkit.Chem import AllChem
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)
salt_remover = SaltRemover.SaltRemover()

def standard_value_to_pic50(val):
    if val <= 0 or pd.isna(val):
        return None
    return -math.log10(val * 1e-9)

def preprocess_molecule(mol):
    if mol is None:
        return None
    frags = rdmolops.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return None
    mol = max(frags, key=lambda m: m.GetNumAtoms())
    mol = salt_remover.StripMol(mol, dontRemoveEverything=True)
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        try:
            mol.UpdatePropertyCache(strict=False)
            Chem.GetSymmSSSR(mol)
        except Exception:
            return None
    return mol

def extract_mmp_like_fragments(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return [], "invalid", smi

    mol = preprocess_molecule(mol)
    if mol is None:
        return [], "preprocessing_failed", smi

    try:
        bonds_to_cut = []
        for bond in mol.GetBonds():
            begin = bond.GetBeginAtom()
            end = bond.GetEndAtom()

            if bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            if bond.IsInRing():
                continue
            if begin.GetAtomicNum() == 1 or end.GetAtomicNum() == 1:
                continue
            if begin.IsInRing() and end.IsInRing():
                continue  # avoid breaking ring-ring bonds

            # Keep only ring-to-chain bonds
            if begin.IsInRing() != end.IsInRing():
                bonds_to_cut.append(bond.GetIdx())

        if not bonds_to_cut:
            return [], "no_valid_cuts", smi

        # Do fragmentation
        fragmented = Chem.FragmentOnBonds(mol, bonds_to_cut, addDummies=True)
        frags = Chem.GetMolFrags(fragmented, asMols=True, sanitizeFrags=True)
        if not frags or len(frags) < 2:
            return [], "too_few_frags", smi

        frag_smis = [(Chem.MolToSmiles(f, isomericSmiles=True), f) for f in frags]
        frag_smis = list(set(frag_smis))  # remove duplicates

        # Identify core (largest frag with ring) and subs (others)
        frag_smis = sorted(frag_smis, key=lambda x: x[1].GetNumAtoms(), reverse=True)
        core = None
        subs = []

        for smi_core, mol_core in frag_smis:
            if mol_core.GetRingInfo().NumRings() > 0:
                core = smi_core
                break

        if core is None:
            return [], "no_ring_core", smi

        for smi_sub, mol_sub in frag_smis:
            if smi_sub == core:
                continue
            subs.append(smi_sub)

        if not subs:
            return [], "no_substituents", smi

        pairs = []
        for sub in subs:
            core_smarts = simplify_fragment_to_smarts(core)
            sub_smarts = simplify_fragment_to_smarts(sub)

            if not core_smarts or not sub_smarts:
                continue

            # Only keep fragments with exactly 1 dummy/attachment point
            # if core_smarts.count('[*]') != 1 or sub_smarts.count('[*]') != 1:
            #     continue

            pairs.append((core_smarts, sub_smarts))

        if not pairs:
            return [], "filtered_all", smi

        return pairs, "has_pairs", smi

    except Exception as e:
        logger.debug(f"Fragmentation error on {smi}: {e}")
        return [], "error", smi

def simplify_fragment_to_smarts(smi: str) -> str:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom.SetIsotope(0)
            atom.SetAtomMapNum(0)
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(0)

    smarts = Chem.MolToSmarts(mol, isomericSmiles=False)
    smarts = smarts.replace('[#0]', '[*]')
    return smarts

@register_task(
    "mmp_rule_learning",
    category="Preprocessing",
    description="Learn BRICS‑based transformation rules from actives"
)
def mmp_rule_learning(config, data=None):
    input_file = config.get("input_file")
    mmp_config = config.get("mmp_rule_learning", {})
    cutoff = mmp_config.get("pic50_cutoff", 7.0)
    filter_rules = mmp_config.get("filter_rules", False)
    min_rule_count = mmp_config.get("min_rule_count", 1)
    min_core_rings = mmp_config.get("min_core_rings", 1)
    min_core_heavy_atoms = mmp_config.get("min_core_heavy_atoms", 6)

    df = pd.read_csv(input_file)
    if "smiles" not in df.columns or "standard_value" not in df.columns:
        raise ValueError("Input CSV must contain 'smiles' and 'standard_value' columns.")

    df["pIC50"] = df["standard_value"].apply(standard_value_to_pic50)
    df = df[df["pIC50"].notna()]
    actives = df[df["pIC50"] >= cutoff]
    logger.info(f"Found {len(actives)} active molecules")

    core_to_subs = defaultdict(set)
    rule_counter = Counter()
    status_counts = Counter()

    with ProcessPoolExecutor() as executor:
        results = list(tqdm(
            executor.map(extract_mmp_like_fragments, actives["smiles"]),
            total=len(actives),
            desc="Learning BRICS rules"
        ))

    total_pairs = 0
    for pairs, status, smi in results:
        status_counts[status] += 1
        for core, sub in pairs:
            try:
                core_mol = Chem.MolFromSmiles(core, sanitize=True)
                if core_mol is None:
                    continue
                if core_mol.GetRingInfo().NumRings() < min_core_rings:
                    continue
                if core_mol.GetNumHeavyAtoms() < min_core_heavy_atoms:
                    continue
            except Exception:
                continue

            core_smarts = simplify_fragment_to_smarts(core)
            sub_smarts = simplify_fragment_to_smarts(sub)

            if core_smarts and sub_smarts:
                core_to_subs[core_smarts].add(sub_smarts)
                total_pairs += 1

    logger.info(f"BRICS extraction status counts: {dict(status_counts)}")
    logger.info(f"Total (core, sub) pairs extracted: {total_pairs}")
    logger.info(f"Unique cores found: {len(core_to_subs)}")

    for core, subs in core_to_subs.items():
        subs = list(subs)
        if len(subs) < 2:
            continue
        for i in range(len(subs)):
            for j in range(i + 1, len(subs)):
                if subs[i].count("[*]") != subs[j].count("[*]"):
                    continue
                rule_counter[(subs[i], subs[j])] += 1
                rule_counter[(subs[j], subs[i])] += 1

    if filter_rules:
        rule_counter = Counter({k: v for k, v in rule_counter.items() if v >= min_rule_count})

    logger.info("Top transformation rules:")
    for (s_from, s_to), cnt in rule_counter.most_common(10):
        logger.info(f"{s_from} => {s_to} | {cnt}")

    scaffold_to_subs = {core: list(subs) for core, subs in core_to_subs.items()}

    out_dir = Path(config.get("output", {}).get("directory", "outputs/mmp_based_transformation"))
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "mmp_rules_brics.json"
    with open(out_file, "w") as f:
        json.dump(scaffold_to_subs, f, indent=2)

    return ("mmp_rules_brics", pd.DataFrame.from_dict(scaffold_to_subs, orient="index"))
