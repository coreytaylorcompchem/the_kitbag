from pipeline.task_registry import register_task
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, rdMolDescriptors
from rdkit.Chem.QED import qed
from pipeline.logger import setup_logger
from pathlib import Path
import pandas as pd

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


def generate_conformer(mol, max_attempts=10):
    mol_h = Chem.AddHs(mol)
    status = AllChem.EmbedMolecule(mol_h, AllChem.ETKDG())
    if status != 0:
        for _ in range(1, max_attempts):
            status = AllChem.EmbedMolecule(mol_h, AllChem.ETKDG())
            if status == 0:
                break
    if status == 0:
        AllChem.UFFOptimizeMolecule(mol_h)
        return mol_h
    return None


@register_task("basic_lipinski_filtering", description="Basic Lipinski Rule of 5 filtering")
def basic_lipinski(config, data=None):
    input_file = config.get("input_file")
    if input_file is None:
        raise ValueError("No input_file specified in config")
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file does not exist: {input_file}")

    df = pd.read_csv(input_file)
    if "smiles" not in df.columns:
        raise ValueError("Input CSV must contain a 'smiles' column.")

    accepted_rows = []

    # Lipinski cutoffs - can override via config.basic_lipinski.{key}_cutoff
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

    # --- New: Canonicalize SMILES ---
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
    # ------------------------------

    logger.debug(f"Number of molecules before filtering: {df.shape[0]}")
    logger.debug(f"Number of molecules after filtering: {filtered_df.shape[0]}")

    logger.debug(f"basic_lipinski_filtering: {len(filtered_df)} molecules passed filters")

    return (input_path.stem, filtered_df)


@register_task("physchem_filtering", description="Physchem filtering with mandatory and optional cutoffs + conformer generation")
def physchem_filtering(config, data=None):
    input_file = config.get("input_file")
    if input_file is None:
        raise ValueError("No input_file specified in config")
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file does not exist: {input_file}")

    df = pd.read_csv(input_file)
    if "smiles" not in df.columns:
        raise ValueError("Input CSV must contain a 'smiles' column.")

    mandatory_cutoffs = config.get("mandatory_cutoffs", {})
    optional_cutoffs = config.get("optional_cutoffs", {})

    accepted_props = []
    mols = []

    for idx, row in df.iterrows():
        smi = row.get("smiles")
        if not smi or pd.isna(smi):
            continue

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            logger.debug(f"Invalid SMILES skipped: {smi}")
            continue

        props = compute_physchem(mol)

        # Check mandatory cutoffs
        if not all(props.get(k, float('inf')) <= mandatory_cutoffs.get(k, float('inf')) for k in MANDATORY_KEYS):
            continue

        # Filter out None from optional cutoffs and apply
        filtered_optional_cutoffs = {k: v for k, v in optional_cutoffs.items() if v is not None}
        if not all(props.get(k, float('inf')) <= v for k, v in filtered_optional_cutoffs.items()):
            continue

        conf_mol = generate_conformer(mol)
        if conf_mol is None:
            continue

        mols.append(conf_mol)
        accepted_props.append(props)

    if not accepted_props:
        logger.debug("No molecules passed filters or conformer generation.")
        empty_df = pd.DataFrame(columns=list(MANDATORY_KEYS) + list(OPTIONAL_KEYS) + ["smiles"])
        return (input_path.stem, empty_df)

    output_df = pd.DataFrame(accepted_props)
    logger.debug(f"physchem_filtering: {len(output_df)} molecules passed filters and conformer generation.")

    # Write output files if specified
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

    return (input_path.stem, output_df)
