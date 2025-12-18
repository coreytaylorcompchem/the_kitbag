import os
import yaml
import csv
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

from backends import get_backend
from modules.protein_preparation import ProteinPreparer
from pipeline.task_registry import get_task
from workflows import register_workflow
from workflows.utils.docking import (
    generate_ligands_csv_from_txt,
    validate_config,
    get_docking_box
)

from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

def run_single_ligand(ligand, backend, config, workflow_steps):
    """
    Run all workflow steps for a single ligand sequentially.
    Returns ligand name for logging.
    """
    for step in workflow_steps:
        task_func = get_task(step)
        if not task_func:
            raise ValueError(f"Workflow step '{step}' is not a registered task.")
        logger.debug(f"→ Running task: {step} for ligand {ligand['name']}")
        task_func(backend, ligand, config)
    return ligand['name']


@register_workflow("induced_fit_docking", description="Perform induced fit docking.")
def run(config_path: str):
    """
    Run the induced-fit docking workflow.
    """
    # ----------------------
    # Load YAML Config
    # ----------------------
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    required_fields = [
        "docking.output_dir",
        "workflow",
        "protein.pdb_path",
        "induced_fit_docking.minimisation.protein_forcefield",
    ]
    validate_config(config, required_fields)

    docking_cfg = config.get("docking", {})
    backend_cfg = config.get("backend", {})

    logger.info("---------------------------------------------")
    logger.info("Loaded YAML configuration:")
    logger.info(f"Backend: {backend_cfg.get('name')} (GPU: {backend_cfg.get('use_gpu')})")
    logger.info(f"Output directory: {docking_cfg.get('output_dir')}")
    logger.info(f"Docking mode: {docking_cfg.get('docking_mode', 'ensemble')}")
    logger.info(f"IFD residue cutoff: {config['induced_fit_docking'].get('residue_distance_cutoff', 5)} Å")
    logger.info(f"Workflow steps: {config.get('workflow', [])}")
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
        raise FileNotFoundError("No ligand input found. Provide 'ligands_txt' or 'ligands_csv'.")

    # ----------------------
    # Backend Initialization
    # ----------------------
    backend_name = backend_cfg["name"]
    backend_kwargs = {k: v for k, v in backend_cfg.items() if k != "name"}
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
    # Docking Box
    # ----------------------
    center, size = get_docking_box(config, protein_preparer)
    config["docking"]["center"] = center
    config["docking"]["size"] = size

    # ----------------------
    # Read Ligands
    # ----------------------
    ligands = []
    with open(ligands_csv_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ligands.append({'name': row['name'], 'smiles': row['smiles']})

    # ----------------------
    # Workflow Steps
    # ----------------------
    workflow_steps = config.get("workflow", [
        "prepare_receptor_pdbqt",
        "standardise_ligand",
        "convert_to_pdbqt",
        "dock",
        "induced_fit_docking"
    ])

    options = config.get("options", {})
    if not options.get("use_conformer_generation", True):
        workflow_steps = [
            step for step in workflow_steps
            if step not in ("generate_conformers", "cluster_conformers", "save_final_conformers")
        ]
        logger.warning("⚠️  Conformer generation steps disabled — skipping conformer-related tasks.")

    # ----------------------
    # Parallel Execution
    # ----------------------
    n_cores = os.cpu_count() or 16
    n_workers = max(1, n_cores - 2)
    logger.info(f"Using {n_workers} parallel workers out of {n_cores} cores.")

    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = [
            executor.submit(run_single_ligand, lig, backend, config, workflow_steps)
            for lig in ligands
        ]
        for future in as_completed(futures):
            try:
                ligand_name = future.result()
                logger.info(f"Completed induced fit docking for ligand: {ligand_name}")
            except Exception as e:
                logger.error(f"❌ Error during IFD for ligand: {e}")

    logger.info("Induced fit docking workflow completed successfully.")
