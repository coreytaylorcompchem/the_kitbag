import yaml
import csv
from pathlib import Path

from backends import get_backend, discover_backends
from modules.protein_preparation import ProteinPreparer
from pipeline.task_registry import get_task
from workflows import register_workflow

from workflows.utils import generate_ligands_csv_from_txt, validate_ligands_csv, get_docking_box, validate_config

from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

# --------------------------
# Main Workflow Function
# --------------------------
@register_workflow("vanilla_docking", description="Prepare, dock and score with no constraints.")
def run(config_path: str):

    # ----------------------
    # Load YAML Config
    # ----------------------

    # Minimum fields required from yaml for successful docking
    required_fields = [
        "docking.output_dir",
        "workflow",
        "docking.final_n_conformers",
        "docking.rmsd_threshold",
        "docking.min_energy_gap",
    ]
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    validate_config(config, required_fields)

    docking_cfg = config.get("docking", {})

    logger.info("---------------------------------------------")
    logger.info("Loaded YAML configuration:")
    logger.info(f"Output directory: {docking_cfg.get('output_dir')}")
    logger.info(f"Conformer generation: max final conformers: {docking_cfg.get('final_n_conformers', 5)}")
    logger.info(f"Conformer generation: RMSD threshold:   {docking_cfg.get('rmsd_threshold', 0.75)}")
    logger.info(f"Conformer generation: Energy gap:       {docking_cfg.get('min_energy_gap', 0.5)}")
    logger.info(f"Docking: output poses: {docking_cfg.get('n_output_binding_modes', 5)}")
    logger.info(f"Docking: exhaustiveness:   {docking_cfg.get('exhaustiveness', 5)}")
    logger.info(f"Workflow steps:   {config.get('workflow', [])}")
    logger.info("---------------------------------------------")

    # ----------------------
    # Output Directory
    # ----------------------
    output_dir = Path.cwd() / docking_cfg['output_dir']
    output_dir.mkdir(parents=True, exist_ok=True)
    config['output_dir'] = output_dir

    # ----------------------
    # Ligand Input Handling
    # ----------------------
    ligands_txt_path = Path(config.get('ligands_txt', 'ligands.txt'))
    ligands_csv_path = Path(config.get('ligands_csv', output_dir / 'ligands.csv'))

    if ligands_txt_path.exists():
        logger.info(f"Found ligands.txt — generating ligands.csv.")
        generate_ligands_csv_from_txt(ligands_txt_path, ligands_csv_path)
    elif ligands_csv_path.exists():
        logger.info(f"Found ligands.csv — using it directly.")
    else:
        raise FileNotFoundError(
            "No ligand input found. Provide either 'ligands_txt' or 'ligands_csv' in the config or directory."
        )
    
    # ----------------------
    # Backend Initialization
    # ----------------------
    # Extract backend details from config
    backend_name = config['backend']['name']
    backend_kwargs = {k: v for k, v in config['backend'].items() if k != 'name'}
    backend = get_backend(backend_name, **backend_kwargs)
    
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
    # Docking box (center, size)
    center, size = get_docking_box(config, protein_preparer)
    config['docking']['center'] = center
    config['docking']['size'] = size

    # ----------------------
    # Ligand CSV Handling
    # ----------------------
    if config.get('generate_ligands_from_list', False):
        if not ligands_txt_path.exists():
            raise FileNotFoundError(f"Ligands SMILES list not found at: {ligands_txt_path}")
        generate_ligands_csv_from_txt(ligands_txt_path, ligands_csv_path)

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
        "standardise_ligand",
        "generate_conformers",
        "cluster_conformers",
        "save_final_conformers",
        "convert_to_pdbqt",
        "dock"
    ])

    # Execute workflow per ligand
    for ligand in ligands:
        logger.info(f"\nProcessing ligand: {ligand['name']}")
        for step in workflow_steps:
            task_func = get_task(step)
            if not task_func:
                raise ValueError(f"Workflow step '{step}' is not a registered task.")
            logger.info(f"Running step: {step}")
            task_func(backend, ligand, config)

    logger.info("\nVanilla docking workflow completed.")
