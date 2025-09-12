import csv
from pathlib import Path
from rdkit import Chem

# --------------------------
# Utility: CSV Generation
# --------------------------
def generate_ligands_csv_from_txt(txt_path: Path, csv_path: Path):
    with open(txt_path, 'r') as f:
        smiles_list = [line.strip() for line in f if line.strip()]

    if not smiles_list:
        raise ValueError(f"No SMILES found in {txt_path}")

    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['name', 'smiles'])
        for i, smi in enumerate(smiles_list, start=1):
            writer.writerow([f"ligand{i}", smi])

    print(f"[INFO] Generated ligands.csv at {csv_path} with {len(smiles_list)} ligands.")

# --------------------------
# Utility: CSV Validation
# --------------------------
def validate_ligands_csv(csv_path: Path):
    if not csv_path.exists():
        raise FileNotFoundError(f"Expected ligands.csv at {csv_path} not found.")
    with open(csv_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        if 'name' not in reader.fieldnames or 'smiles' not in reader.fieldnames:
            raise ValueError(f"ligands.csv must have 'name' and 'smiles' columns.")

# --------------------------
# Utility: Docking Center
# --------------------------
def get_docking_box(config, protein_preparer):
    """
    Determine the docking box center and size from config.
    """
    docking_cfg = config.get('docking', {})

    # -------- Get docking center --------
    if 'center' in docking_cfg:
        center = tuple(docking_cfg['center'])

    elif 'ref_ligand_path' in docking_cfg:
        ref_path = docking_cfg['ref_ligand_path']
        ligand_mol = Chem.MolFromMolFile(ref_path)
        if ligand_mol is None:
            raise ValueError(f"[ERROR] Reference ligand file invalid: {ref_path}")
        conf = ligand_mol.GetConformer()
        pts = [list(conf.GetAtomPosition(i)) for i in range(ligand_mol.GetNumAtoms())]
        center = tuple(sum(x) / len(x) for x in zip(*pts))

    elif docking_cfg.get('use_crystal_ligand', False):
        ref_lig_path = protein_preparer.reference_ligand_path
        ligand_mol = Chem.MolFromPDBFile(str(ref_lig_path))
        if ligand_mol is None:
            raise ValueError("[ERROR] Failed to read crystallized ligand from protein.")
        conf = ligand_mol.GetConformer()
        pts = [list(conf.GetAtomPosition(i)) for i in range(ligand_mol.GetNumAtoms())]
        center = tuple(sum(x) / len(x) for x in zip(*pts))

    else:
        raise ValueError(
            "[ERROR] Docking center not specified. Provide one of:\n"
            "  - center: [x, y, z]\n"
            "  - ref_ligand_path: /path/to/ligand\n"
            "  - use_crystal_ligand: true"
        )

    # -------- Get docking box size --------
    if "size" not in docking_cfg:
        raise ValueError("[ERROR] Docking 'size' must be specified in config['docking'].")

    size = tuple(docking_cfg["size"])

    return center, size

# --------------------------
# Utility: Validate yaml 
# --------------------------

def validate_config(config, required_fields):
    missing_fields = []
    for field in required_fields:
        parts = field.split('.')
        value = config
        for part in parts:
            if isinstance(value, dict) and part in value:
                value = value[part]
            else:
                value = None
                break
        if value in (None, [], ''):
            missing_fields.append(field)
    
    if missing_fields:
        print(f"[WARNING] Required config fields are missing: {missing_fields}")
