from pipeline.task_registry import register_task
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, rdMolDescriptors
from rdkit.Chem.QED import qed
from pipeline.logger import setup_logger

from tqdm import tqdm
from pathlib import Path
import pandas as pd
from collections import Counter

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

MANDATORY_KEYS = ["mw", "hbd", "hba"]
OPTIONAL_KEYS = {"logp", "rotatable_bonds", "tpsa", "qed", "stereocenters"}


def compute_physchem(mol):
    return {
        "mol": mol,
        "smiles": Chem.MolToSmiles(mol),
        "mw": Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "hbd": Descriptors.NumHDonors(mol),
        "hba": Descriptors.NumHAcceptors(mol),
        "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
        "tpsa": Descriptors.TPSA(mol),
        "qed": qed(mol),
        "stereocenters": rdMolDescriptors.CalcNumAtomStereoCenters(mol),
    }


def apply_mandatory_filters(props, cutoffs):
    return all(props[k] <= cutoffs[k] for k in MANDATORY_KEYS)


def apply_optional_filters(props, cutoffs):
    for key, threshold in cutoffs.items():
        val = props.get(key)
        if val is None or val > threshold:
            return False
    return True

from rdkit.Chem import rdDistGeom

def generate_conformer(mol, max_attempts=10):
    mol_h = Chem.AddHs(mol)

    status = AllChem.EmbedMolecule(mol_h, AllChem.ETKDG())
    if status != 0:
        return None, "embed_failed"

    AllChem.UFFOptimizeMolecule(mol_h)
    return mol_h, None


@register_task("basic_lipinski_filtering", category="Filtering", description="Basic Lipinski Rule of 5 filtering")
def basic_lipinski(config, data=None):
    input_file = config.get("input_file")
    if input_file is None:
        raise ValueError("No input_file specified in config")
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file does not exist: {input_file}")

    df = pd.read_csv(input_file)
    if data is not None and "df" in data:
        df = data["df"].copy()
    else:
        df = pd.read_csv(input_file)
    if "smiles" not in df.columns:
        raise ValueError("Input CSV must contain a 'smiles' column.")

    accepted_rows = []

    # Lipinski cutoffs
    mw_cutoff = config.get("basic_lipinski", {}).get("mw_cutoff", 500)
    hbd_cutoff = config.get("basic_lipinski", {}).get("hbd_cutoff", 5)
    hba_cutoff = config.get("basic_lipinski", {}).get("hba_cutoff", 10)
    logp_cutoff = config.get("basic_lipinski", {}).get("logp_cutoff", 5.0)

    for idx, row in df.iterrows():
        smi = row.get("smiles")
        if not smi or pd.isna(smi):
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            logger.debug(f"Failed to parse SMILES: {smi}")
            continue

        mw = Descriptors.MolWt(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        logp = Descriptors.MolLogP(mol)

        if mw <= mw_cutoff and hbd <= hbd_cutoff and hba <= hba_cutoff and logp <= logp_cutoff:
            accepted_rows.append(row.to_dict())

    if not accepted_rows:
        logger.debug("No molecules in this batch passed filters/conformer generation.")
        empty_df = pd.DataFrame(columns=["smiles"] + MANDATORY_KEYS + list(OPTIONAL_KEYS))
        return (Path(input_file).stem, empty_df)

    filtered_df = pd.DataFrame(accepted_rows)

    def canonicalize_smi(smi):
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                return Chem.MolToSmiles(mol, canonical=True)
        except Exception:
            pass
        return None

    filtered_df["smiles"] = filtered_df["smiles"].astype(str).str.strip()
    filtered_df["smiles"] = filtered_df["smiles"].apply(canonicalize_smi)
    filtered_df = filtered_df[filtered_df["smiles"].notna()].reset_index(drop=True)

    logger.debug(f"Number of molecules before filtering: {df.shape[0]}")
    logger.debug(f"Number of molecules after filtering: {filtered_df.shape[0]}")

    logger.debug(f"basic_lipinski_filtering: {len(filtered_df)} molecules passed filters")

    return (input_path.stem, filtered_df)

@register_task("physchem_filtering", category="Filtering", description="Physchem filtering with mandatory and optional cutoffs + conformer generation")
def physchem_filtering(config, data=None):
    input_file = config.get("input_file")
    if input_file is None:
        raise ValueError("No input_file specified in config")
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file does not exist: {input_file}")

    if data is not None and "df" in data:
        df = data["df"].copy()
    else:
        df = pd.read_csv(input_file)

    if "smiles" not in df.columns:
        raise ValueError("Input CSV must contain a 'smiles' column.")

    mandatory_cutoffs = config.get("mandatory_cutoffs", {})
    optional_cutoffs = config.get("optional_cutoffs", {})

    accepted_rows = []
    mols = []

    use_progress = "chunk_size" not in config 

    iterator = df.iterrows()
    if use_progress:
        iterator = tqdm(iterator, total=len(df), desc="Processing molecules")

    for idx, row in iterator:
        smi = row.get("smiles")
        if not smi or pd.isna(smi):
            continue

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            logger.debug(f"Invalid SMILES skipped: {smi}")
            continue

        props = compute_physchem(mol)

        if not all(props.get(k, float("inf")) <= mandatory_cutoffs.get(k, float("inf")) for k in MANDATORY_KEYS):
            continue

        filtered_optional_cutoffs = {k: v for k, v in optional_cutoffs.items() if v is not None}
        if not all(props.get(k, float("inf")) <= v for k, v in filtered_optional_cutoffs.items()):
            continue

        conf_mol, conf_fail = generate_conformer(mol)
        
        fail_counts = Counter()

        if conf_mol is None:
            fail_counts[f"conformer_{conf_fail}"] += 1

            # optional: log examples but don’t spam
            if fail_counts[f"conformer_{conf_fail}"] <= 3:
                logger.debug(
                    f"Conformer failed ({conf_fail}) for SMILES: {smi}"
                )
            continue


        row_data = row.to_dict()
        combined = {**row_data, **props}

        accepted_rows.append(combined)
        mols.append(conf_mol)

        if use_progress:
            iterator.set_postfix({"accepted": len(accepted_rows)})

    if not accepted_rows:
        logger.debug("No molecules passed filters or conformer generation.")
        all_columns = list(df.columns) + list(MANDATORY_KEYS) + list(OPTIONAL_KEYS)
        all_columns = list(dict.fromkeys(all_columns))  # remove duplicates
        empty_df = pd.DataFrame(columns=all_columns)
        return (input_path.stem, empty_df)

    output_df = pd.DataFrame(accepted_rows)
    logger.debug(f"physchem_filtering: {len(output_df)} molecules passed filters and conformer generation.")

    output_dir = Path(config.get("output", {}).get("directory", "."))
    output_dir.mkdir(parents=True, exist_ok=True)

    base_filename = config.get("output", {}).get("filename", input_path.stem)

    smi_file = output_dir / f"{base_filename}.smi"
    sdf_file = output_dir / f"{base_filename}.sdf"

    try:
        with open(smi_file, "w") as f:
            for smi in output_df["smiles"].dropna():
                f.write(smi + "\n")
        logger.debug(f"SMILES written to {smi_file}")
    except Exception as e:
        logger.error(f"Failed to write SMILES file {smi_file}: {e}")

    try:
        writer = Chem.SDWriter(str(sdf_file))
        for mol in mols:
            writer.write(mol)
        writer.close()
        logger.debug(f"SDF file written to {sdf_file}")
    except Exception as e:
        logger.error(f"Failed to write SDF file {sdf_file}: {e}")

    total = len(df)
    passed = len(output_df)
    failed = total - passed

    logger.info(
        f"physchem_filtering summary: "
        f"total={total}, passed={passed}, failed={failed}"
    )

    if fail_counts:
        for reason, count in fail_counts.most_common():
            logger.info(f"  {reason}: {count}")

    return (input_path.stem, output_df)