import yaml
import csv
from pathlib import Path
from rdkit import Chem

from backends.gnina import GninaBackend
from modules.protein_preparation import ProteinPreparer
from docking_task_registry import get_task
from workflows import register_workflow

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
def get_docking_center(config, protein_preparer):
    docking_cfg = config['docking']

    if 'center' in docking_cfg:
        return tuple(docking_cfg['center'])

    elif 'ref_ligand_path' in docking_cfg:
        ligand_mol = Chem.MolFromMolFile(docking_cfg['ref_ligand_path'])
        if ligand_mol is None:
            raise ValueError(f"Reference ligand file invalid: {docking_cfg['ref_ligand_path']}")
        conf = ligand_mol.GetConformer()
        pts = [list(conf.GetAtomPosition(i)) for i in range(ligand_mol.GetNumAtoms())]
        return tuple(sum(x)/len(x) for x in zip(*pts))

    elif docking_cfg.get('use_crystal_ligand', False):
        ref_lig_path = protein_preparer.reference_ligand_path
        ligand_mol = Chem.MolFromPDBFile(str(ref_lig_path))
        if ligand_mol is None:
            raise ValueError("Failed to read crystallized ligand from protein.")
        conf = ligand_mol.GetConformer()
        pts = [list(conf.GetAtomPosition(i)) for i in range(ligand_mol.GetNumAtoms())]
        return tuple(sum(x)/len(x) for x in zip(*pts))

    else:
        raise ValueError("Docking center not specified. Provide `center`, `ref_ligand_path`, or `use_crystal_ligand: true`.")

# --------------------------
# Main Workflow Function
# --------------------------
@register_workflow("vanilla_docking", description="Preparation and docking using Gnina.")
def run(config_path: str):
    # ----------------------
    # Load YAML Config
    # ----------------------
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    print("[INFO] Loaded YAML configuration:")
    print(f"  - Output directory: {config.get('output_dir')}")
    print(f"  - Final conformers: {config.get('final_n_conformers', 5)}")
    print(f"  - RMSD threshold:   {config.get('rmsd_threshold', 0.75)}")
    print(f"  - Energy gap:       {config.get('min_energy_gap', 0.5)}")
    print(f"  - Workflow steps:   {config.get('workflow', [])}")

    # ----------------------
    # Output Directory
    # ----------------------
    output_dir = Path.cwd() / config['output_dir']
    output_dir.mkdir(parents=True, exist_ok=True)

    # ----------------------
    # Ligand Input Handling
    # ----------------------
    ligands_txt_path = Path(config.get('ligands_txt', 'ligands.txt'))
    ligands_csv_path = Path(config.get('ligands_csv', output_dir / 'ligands.csv'))

    if ligands_txt_path.exists():
        print(f"[INFO] Found ligands.txt — generating ligands.csv.")
        generate_ligands_csv_from_txt(ligands_txt_path, ligands_csv_path)

    elif ligands_csv_path.exists():
        print(f"[INFO] Found ligands.csv — using it directly.")

    else:
        raise FileNotFoundError("No ligand input found. Provide either 'ligands_txt' or 'ligands_csv' in the config or directory.")
    
    # ----------------------
    # Backend Initialization
    # ----------------------
    backend = GninaBackend(
        gnina_executable=config['backend']['gnina_path'],
        use_gpu=config['backend'].get('use_gpu', True)
    )
    
    # ----------------------
    # Protein Preparation
    # ----------------------
    protein_preparer = ProteinPreparer(
        pdb_path=Path(config['protein']['pdb_path']),
        work_dir=output_dir,
        pH=config['protein'].get('pH', 7.4)
    )
    receptor_pdbqt = protein_preparer.prepare()
    backend.cache["receptor_pdbqt"] = receptor_pdbqt

    # ----------------------
    # Docking Center
    # ----------------------
    center = get_docking_center(config, protein_preparer)
    size = tuple(config['docking']['size'])

    # Inject back into config for tasks
    config['docking']['center'] = center
    config['docking']['size'] = size

    # ----------------------
    # Ligand CSV Handling
    # ----------------------
    if config.get('generate_ligands_from_list', False):
        ligands_txt_path = Path(config.get('ligands_txt', 'ligands.txt'))

        if not ligands_txt_path.exists():
            raise FileNotFoundError(f"Ligands SMILES list not found at: {ligands_txt_path}")
        
        # Default output CSV path if not specified
        ligands_csv_path = Path(config.get('ligands_csv', 'ligands.csv'))
        generate_ligands_csv_from_txt(ligands_txt_path, ligands_csv_path)

    else:
        ligands_csv_path = Path(config.get('ligands_csv', 'ligands.csv'))
        if not ligands_csv_path.exists():
            raise FileNotFoundError(f"Ligands CSV not found at: {ligands_csv_path}")

    # ----------------------
    # Read Ligands
    # ----------------------
    ligands = []
    with open(ligands_csv_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ligands.append({'name': row['name'], 'smiles': row['smiles']})

    # ----------------------
    # Execute Workflow Steps
    # ----------------------
    workflow_steps = config.get("workflow", [
        "standardize_ligand",
        "generate_conformers",
        "cluster_conformers",
        "save_final_conformers",
        "convert_to_pdbqt",
        "dock"
    ])

    for ligand in ligands:
        print(f"\n[INFO] Processing ligand: {ligand['name']}")
        for step in workflow_steps:
            task_func = get_task(step)
            if not task_func:
                raise ValueError(f"Workflow step '{step}' is not a registered task.")
            print(f"[INFO] Running step: {step}")
            task_func(backend, ligand, config)

    print("\n[INFO] Vanilla docking workflow completed.")
