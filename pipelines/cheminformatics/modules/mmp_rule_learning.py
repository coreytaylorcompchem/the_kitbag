from pipeline.task_registry import register_task
import math, json
from pathlib import Path
import pandas as pd
from tqdm import tqdm
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor

from rdkit import Chem
from rdkit.Chem import BRICS, rdmolops
from rdkit.Chem import SaltRemover

from pipeline.logger import setup_logger
logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def standard_value_to_pic50(val):
    if val <= 0 or pd.isna(val):
        return None
    return -math.log10(val * 1e-9)

# Salt remover
salt_remover = SaltRemover.SaltRemover()

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

def normalize_dummy_atoms(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return smi  # fallback if parsing fails

    rw = Chem.RWMol(mol)
    dummy_atoms = sorted([a for a in rw.GetAtoms() if a.GetAtomicNum() == 0], key=lambda a: a.GetIdx())
    # Set isotope number to 1, 2, 3... as dummy atom labels
    for i, atom in enumerate(dummy_atoms, start=1):
        atom.SetIsotope(i)
    # Generate canonical smiles preserving isotopes (dummy labels)
    smi_norm = Chem.MolToSmiles(rw, canonical=True)
    return smi_norm

def extract_brics_fragments(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return [], "invalid", smi
    mol = preprocess_molecule(mol)
    if mol is None:
        return [], "preprocessing_failed", smi

    try:
        frag_smis = list(BRICS.BRICSDecompose(mol))
        if not frag_smis:
            return [], "no_brics_fragments", smi

        # Convert to Mol objects and count atoms
        frags = []
        for f_smi in frag_smis:
            f_mol = Chem.MolFromSmiles(f_smi)
            if f_mol is None:
                continue
            frags.append((f_smi, f_mol.GetNumAtoms()))

        if len(frags) < 2:
            return [], "too_few_fragments", smi

        # Pick core = largest fragment, others are substituents
        frags = sorted(frags, key=lambda x: x[1], reverse=True)
        core_smi, _ = frags[0]
        subs = [f[0] for f in frags[1:]]

        pairs = [(core_smi, sub) for sub in subs if sub != core_smi]
        return list(set(pairs)), "has_brics_pairs", smi

    except Exception as e:
        logger.debug(f"BRICS extraction error on {smi}: {e}")
        return [], "error", smi

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
            executor.map(extract_brics_fragments, actives["smiles"]),
            total=len(actives),
            desc="Learning BRICS rules"
        ))

    total_pairs = 0
    for pairs, status, smi in results:
        status_counts[status] += 1
        for core, sub in pairs:
            # Normalize dummy atoms here before storing
            core_norm = normalize_dummy_atoms(core)
            sub_norm = normalize_dummy_atoms(sub)
            core_to_subs[core_norm].add(sub_norm)
            total_pairs += 1

    logger.info(f"BRICS extraction status counts: {dict(status_counts)}")
    logger.info(f"Total (core, sub) pairs extracted: {total_pairs}")
    logger.info(f"Unique cores found: {len(core_to_subs)}")

    # Derive transformation rules (core → substituent pairs)
    for core, subs in core_to_subs.items():
        subs = list(subs)
        if len(subs) < 2:
            continue
        for i in range(len(subs)):
            for j in range(i + 1, len(subs)):
                rule_counter[(subs[i], subs[j])] += 1
                rule_counter[(subs[j], subs[i])] += 1

    if filter_rules:
        rule_counter = Counter({k: v for k, v in rule_counter.items() if v >= min_rule_count})

    logger.info("Top transformation rules:")
    for (s_from, s_to), cnt in rule_counter.most_common(10):
        logger.info(f"{s_from} => {s_to} | {cnt}")

    # Save rule mapping: core → list of substituents
    scaffold_to_subs = {core: list(subs) for core, subs in core_to_subs.items()}

    out_dir = Path(config.get("output", {}).get("directory", "outputs/mmp_based_transformation"))
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "mmp_rules_brics.json"
    with open(out_file, "w") as f:
        json.dump(scaffold_to_subs, f, indent=2)

    return ("mmp_rules_brics", pd.DataFrame.from_dict(scaffold_to_subs, orient="index"))
