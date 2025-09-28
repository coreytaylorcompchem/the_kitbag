from pipeline.task_registry import register_task
import gc
import math
import json
from pathlib import Path
import pandas as pd
from tqdm import tqdm
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor

from rdkit import Chem
from rdkit.Chem import rdMMPA, AllChem, rdmolops, SaltRemover

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)


def standard_value_to_pic50(val):
    if val <= 0 or pd.isna(val):
        return None
    return -math.log10(val * 1e-9)  # Convert nM to pIC50


# Initialize salt remover once (common salts)
salt_remover = SaltRemover.SaltRemover()

def preprocess_molecule(mol):
    if mol is None:
        return None
    
    # Remove salts/fragments, keep largest fragment
    frags = rdmolops.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return None
    mol = max(frags, key=lambda m: m.GetNumAtoms())
    
    # Remove salts using SaltRemover
    mol = salt_remover.StripMol(mol, dontRemoveEverything=True)
    
    # Try sanitizing mol robustly
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        try:
            mol.UpdatePropertyCache(strict=False)
            Chem.GetSymmSSSR(mol)
        except Exception:
            return None
    
    return mol

def extract_mmp_pairs(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return [], "invalid", smi

    # Preprocess molecule
    mol = preprocess_molecule(mol)
    if mol is None:
        return [], "preprocessing_failed", smi

    try:
        # Try fragmentation with maxCuts=2 first
        cuts = rdMMPA.FragmentMol(mol, maxCuts=2, maxCutBonds=60, pattern="[*:1]~[*:2]")
        if not cuts:
            # Fallback: try maxCuts=1 if no cuts found
            cuts = rdMMPA.FragmentMol(mol, maxCuts=1, maxCutBonds=60, pattern="[*:1]~[*:2]")
        if not cuts:
            return [], "no_mmp_cuts", smi

        pairs = []
        for core_mol, side_chain in cuts:
            if core_mol is None or side_chain is None:
                continue
            core_smiles = Chem.MolToSmiles(core_mol, canonical=True)
            sub_smiles = Chem.MolToSmiles(side_chain, canonical=True)
            pairs.append((core_smiles, sub_smiles))

        if not pairs:
            return [], "no_mmp_pairs", smi

        return list(set(pairs)), "has_mmp_pairs", smi

    except Exception as e:
        logger.debug(f"MMP extraction error on {smi}: {e}")
        return [], "error", smi



@register_task(
    "mmp_rule_learning",
    category="Preprocessing",
    description="Learn MMPA-based transformation rules from actives"
)
def mmp_rule_learning(config, data=None):
    input_file = config.get("input_file")
    mmp_config = config.get("mmp_rule_learning", {})
    cutoff = mmp_config.get("pic50_cutoff", 7.0)
    filter_rules = mmp_config.get("filter_rules", False)
    min_rule_count = mmp_config.get("min_rule_count", 2)

    if not Path(input_file).exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

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
    uncuttable_smiles = []

    with ProcessPoolExecutor() as executor:
        results = list(tqdm(
            executor.map(extract_mmp_pairs, actives["smiles"]),
            total=len(actives),
            desc="Learning MMP rules"
        ))

    total_pairs = 0
    for pairs, status, smi in results:
        status_counts[status] += 1
        if status != "has_mmp_pairs":
            uncuttable_smiles.append((smi, status))
        total_pairs += len(pairs)
        for core, sub in pairs:
            core_to_subs[core].add(sub)

    logger.debug(f"MMP extraction stats: {dict(status_counts)}")
    logger.debug(f"Total (core, sub) pairs extracted: {total_pairs}")
    logger.debug(f"Unique cores found: {len(core_to_subs)}")

    # Filter cores with >1 substituent
    cores_with_multiple_subs = [(core, subs) for core, subs in core_to_subs.items() if len(subs) > 1]
    logger.debug(f"Cores with >1 substituent (can generate rules): {len(cores_with_multiple_subs)}")

    sample_out = Path("outputs/mmp_based_transformation/core_sub_examples.txt")
    sample_out.parent.mkdir(parents=True, exist_ok=True)
    with open(sample_out, "w") as f:
        for core, subs in cores_with_multiple_subs[:20]:
            f.write(f"Core: {core}\n")
            for sub in subs:
                f.write(f"  Sub: {sub}\n")
            f.write("\n")
    logger.info(f"Sample core/substituent examples saved to: {sample_out}")

    # Log SMILES that could not be cut
    fail_out = Path("outputs/mmp_based_transformation/uncuttable_smiles.txt")
    with open(fail_out, "w") as f:
        for smi, status in uncuttable_smiles:
            f.write(f"{smi}  # {status}\n")
    logger.debug(f"Uncuttable SMILES saved to: {fail_out}")

    # Derive transformation rules
    rule_count = 0
    for core, subs in core_to_subs.items():
        subs = list(subs)
        if len(subs) < 2:
            continue
        for i in range(len(subs)):
            for j in range(i + 1, len(subs)):
                rule_counter[(subs[i], subs[j])] += 1
                rule_counter[(subs[j], subs[i])] += 1
                rule_count += 2

    logger.debug(f"Total transformation rules generated: {rule_count}")
    logger.debug(f"Total unique rules: {len(rule_counter)}")

    if filter_rules:
        before = len(rule_counter)
        rule_counter = Counter({k: v for k, v in rule_counter.items() if v >= min_rule_count})
        after = len(rule_counter)
        logger.info(f"Filtered rules: {before} → {after} (min_rule_count={min_rule_count})")

    logger.info(f"Top 10 most frequent rules:")
    for (sub_from, sub_to), count in rule_counter.most_common(10):
        logger.info(f"{sub_from} => {sub_to}  |  count: {count}")

    # Save rules to JSON
    scaffold_to_subs = {core: list(subs) for core, subs in core_to_subs.items()}

    out_dir = Path(config.get("output", {}).get("directory", "outputs/mmp_based_transformation"))
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "mmp_rules.json"

    with open(out_file, "w") as f:
        json.dump(scaffold_to_subs, f, indent=2)

    logger.info(f"Transformation rules saved to: {out_file}")

    return ("mmp_rules", pd.DataFrame.from_dict(scaffold_to_subs, orient="index"))
