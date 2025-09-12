import yaml
import csv
from pathlib import Path

from backends.gnina import GninaBackend
from modules.protein_preparation import ProteinPreparer
from docking_task_registry import get_task
from workflows import register_workflow

from workflows.utils import generate_ligands_csv_from_txt, validate_ligands_csv, get_docking_box, validate_config

# --------------------------
# Main Workflow Function
# --------------------------
@register_workflow("vanilla_docking", description="Preparation and docking using Gnina.")
def run(config_path: str):

    # ----------------------
    # Load YAML Config
    # ----------------------

    # Minimum fields required from yaml for successful docking
    required_fields = [
        "output_dir",
        "workflow",
        "docking.final_n_conformers",
        "docking.rmsd_threshold",
        "docking.min_energy_gap",
    ]
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    validate_config(config, required_fields) # check to make sure all minimum steps are in the workflow

    docking_cfg = config.get("docking", {})

    print("[INFO] Loaded YAML configuration:")
    print(f"  - Output directory: {config.get('output_dir')}")
    print(f"  - Final conformers: {docking_cfg.get('final_n_conformers', 5)}")
    print(f"  - RMSD threshold:   {docking_cfg.get('rmsd_threshold', 0.75)}")
    print(f"  - Energy gap:       {docking_cfg.get('min_energy_gap', 0.5)}")
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
    center, size = get_docking_box(config, protein_preparer)

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
