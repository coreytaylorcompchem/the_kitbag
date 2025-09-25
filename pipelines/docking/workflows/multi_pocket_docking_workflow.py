import yaml
import csv
from pathlib import Path

from backends import get_backend
from modules.protein_preparation import ProteinPreparer
from pipeline.task_registry import get_task
from workflows import register_workflow

from workflows.utils.docking import (
    generate_ligands_csv_from_txt,
    summarise_docking_results,
    validate_config,
)
from workflows.utils.fpocket import (
    run_fpocket,
    extract_pocket_centers,
    plot_multi_pocket_scores
)

from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)


@register_workflow("multi_pocket_docking", description="Dock ligands into top N pockets detected by fpocket.")
def run(config_path: str):
    # ----------------------
    # Load YAML Config
    # ----------------------
    required_fields = [
        "docking.output_dir",
        "workflow",
        "docking.final_n_conformers",
        "docking.rmsd_threshold",
        "docking.min_energy_gap",
        "docking.n_pockets_to_use"
    ]

    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    validate_config(config, required_fields)

    docking_cfg = config["docking"]
    output_dir = Path(docking_cfg["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    config["output_dir"] = output_dir

    logger.info("---------------------------------------------")
    logger.info("Loaded YAML configuration:")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Docking: using top {docking_cfg['n_pockets_to_use']} pockets from fpocket.")
    logger.info(f"Workflow steps: {config.get('workflow', [])}")
    logger.info("---------------------------------------------")

    # ----------------------
    # Ligand Input
    # ----------------------
    ligands_txt_path = Path(config.get('ligands_txt', 'ligands.txt'))
    ligands_csv_path = Path(config.get('ligands_csv', output_dir / 'ligands.csv'))

    if ligands_txt_path.exists():
        logger.info(f"Found ligands.txt — generating ligands.csv.")
        generate_ligands_csv_from_txt(ligands_txt_path, ligands_csv_path)
    elif ligands_csv_path.exists():
        logger.info(f"Found ligands.csv — using it directly.")
    else:
        raise FileNotFoundError("No ligand input found.")

    # ----------------------
    # Backend
    # ----------------------
    backend_name = config['backend']['name']
    backend_kwargs = {k: v for k, v in config['backend'].items() if k != 'name'}
    backend = get_backend(backend_name, **backend_kwargs)

    # ----------------------
    # Protein Prep
    # ----------------------
    protein_preparer = ProteinPreparer(
        pdb_path=Path(config['protein']['pdb_path']),
        work_dir=output_dir,
        pH=config['protein'].get('pH', 7.4)
    )
    receptor_pdbqt = protein_preparer.prepare()
    backend.cache["receptor_pdbqt"] = receptor_pdbqt

    # ----------------------
    # Run fpocket and extract pockets
    # ----------------------
    fpocket_output_dir = run_fpocket(protein_preparer.protonated_pdb, output_dir)
    pocket_centers = extract_pocket_centers(fpocket_output_dir, top_n=docking_cfg["n_pockets_to_use"])

    # ----------------------
    # Read Ligands
    # ----------------------
    ligands = []
    with open(ligands_csv_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ligands.append({'name': row['name'], 'smiles': row['smiles']})

    # ----------------------
    # Multi-pocket Docking
    # ----------------------
    for i, center in enumerate(pocket_centers):
        logger.info(f"\n--- Docking into pocket {i+1}/{len(pocket_centers)} ---")
        config['docking']['center'] = center
        config['docking']['size'] = [20, 20, 20]  # Could be dynamic later

        for ligand in ligands:
            logger.info(f"Processing ligand: {ligand['name']} for pocket {i+1}")
            for step in config.get("workflow", []):
                task_func = get_task(step)
                if not task_func:
                    raise ValueError(f"Workflow step '{step}' is not a registered task.")
                task_func(backend, ligand, config, pocket_id=i+1)  # Pass pocket_id if needed

    # ----------------------
    # Summarize & Plot
    # ----------------------
    results_csv = output_dir / "multi_pocket_docking_summary.csv"
    summarise_docking_results(output_dir, results_csv)
    if config.get("plot_docking_results", True):
        plot_multi_pocket_scores(results_csv, output_dir)
    else:
        logger.info("Skipping result plotting.")

    logger.info("Multi-pocket docking workflow completed.")