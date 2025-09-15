from pipeline.task_registry import register_task
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, rdMolDescriptors
from rdkit.Chem.QED import qed
from pipeline.logger import setup_logger
from pathlib import Path
import pandas as pd

logger = setup_logger(__name__)

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
    smiles_list = df["smiles"].dropna().unique().tolist()
    logger.debug(f"[physchem_filtering] Loaded {len(smiles_list)} SMILES")

    # Extract cutoff configs - keys and default values
    mandatory_cutoffs = {
        key: config.get("basic_lipinski", {}).get(f"{key}_cutoff", {
            "mw": 500,
            "logp": 5.0,
            "hbd": 5,
            "hba": 10,
        }[key])
        for key in MANDATORY_KEYS
    }

    optional_cutoffs = {
        key: config.get("physchem", {}).get(f"{key}_cutoff")
        for key in OPTIONAL_KEYS
        if config.get("physchem", {}).get(f"{key}_cutoff") is not None
    }

    accepted = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        props = compute_physchem(mol)
        if not apply_mandatory_filters(props, mandatory_cutoffs):
            continue
        if optional_cutoffs and not apply_optional_filters(props, optional_cutoffs):
            continue
        confmol = generate_conformer(mol)
        if confmol is None:
            continue
        props["mol"] = confmol
        accepted.append(props)

    if not accepted:
        logger.warning("No molecules in this batch passed filters/conformer generation.")
        return {"mols": [], "df": pd.DataFrame()}

    # Save SMILES and SDF
    output_dir = Path(config.get("output", {}).get("directory", "outputs/physchem"))
    output_dir.mkdir(parents=True, exist_ok=True)

    smiles_out = output_dir / "filtered_mols.smi"
    sdf_out = output_dir / "filtered_mols.sdf"

    with open(smiles_out, "w") as f_sm:
        for p in accepted:
            f_sm.write(p["smiles"] + "\n")

    writer = Chem.SDWriter(str(sdf_out))
    for p in accepted:
        writer.write(p["mol"])
    writer.close()

    logger.debug(f"Saved {len(accepted)} molecules to {smiles_out} and {sdf_out}")

    df_rows = []
    for p in accepted:
        row = {k: p.get(k) for k in ["smiles"] + MANDATORY_KEYS + list(OPTIONAL_KEYS)}
        df_rows.append(row)

    return {
        "mols": [p["mol"] for p in accepted],
        "df": pd.DataFrame(df_rows)
    }
