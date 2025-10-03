import yaml
import csv
from pathlib import Path
import pandas as pd

from backends import get_backend
from modules.protein_preparation import ProteinPreparer
from pipeline.task_registry import get_task
from workflows import register_workflow

from workflows.utils.docking import (
    generate_ligands_csv_from_txt,
    summarise_docking_results,
    validate_config,
    extract_scores_from_docking_output
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
    # Handle conformer generation flag
    # ----------------------
    options = config.get("options", {})
    use_conformer_generation = options.get("use_conformer_generation", True)

    workflow_steps = config.get("workflow", [])
    if not use_conformer_generation:
        logger.warning("⚠️  Conformer generation steps disabled - skipping.")
        workflow_steps = [
            step for step in workflow_steps
            if step not in ("generate_conformers", "cluster_conformers", "save_final_conformers")
        ]

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
        logger.info(f">>>>>>>>>> Docking into pocket {i+1}/{len(pocket_centers)} <<<<<<<<<<")
        config['docking']['center'] = center
        config['docking']['size'] = [20, 20, 20]  # TODO: make dynamic later

        for ligand in ligands:
            logger.info(f"Processing ligand: {ligand['name']} for pocket {i+1}")
            docking_outputs = []

            for step in workflow_steps:
                task_func = get_task(step)
                if not task_func:
                    raise ValueError(f"Workflow step '{step}' is not a registered task.")
                result = task_func(backend, ligand, config, pocket_id=i+1)

                if isinstance(result, list):
                    docking_outputs.extend(result)
                elif result:
                    docking_outputs.append(result)

            logger.debug(f"docking_outputs = {docking_outputs}")

            if docking_outputs:
                pocket_dir = output_dir / f"pocket_{i+1}"
                pocket_dir.mkdir(parents=True, exist_ok=True)

                score_rows = []

                for entry in docking_outputs:
                    pdbqt = entry["pdbqt"]
                    conformer = entry["conformer"]
                    docked_outputs = entry["docked_output"]

                    if not isinstance(docked_outputs, (list, tuple)):
                        docked_outputs = [docked_outputs]

                    for sdf_path in docked_outputs:
                        if isinstance(sdf_path, list):
                            if len(sdf_path) == 0:
                                logger.warning(f"No docking outputs found for ligand {ligand['name']}")
                                continue
                            sdf_path = sdf_path[0]

                        score_data = extract_scores_from_docking_output(sdf_path)

                        for score_entry in score_data:
                            score_rows.append({
                                "ligand": ligand["name"],
                                "conformer": conformer,
                                "pose_idx": score_entry["pose_idx"],
                                "score": score_entry["score"],
                                "pose_rank": score_entry["pose_rank"],
                                "pdbqt": pdbqt,
                                "docked_sdf": str(sdf_path),
                                "pocket": i + 1  # use 1-based indexing
                            })

                # Save scores for this pocket
                score_csv_path = pocket_dir / "docking_scores.csv"
                if score_csv_path.exists() and score_csv_path.stat().st_size > 0:
                    existing_df = pd.read_csv(score_csv_path)
                    df = pd.concat([existing_df, pd.DataFrame(score_rows)], ignore_index=True)
                else:
                    df = pd.DataFrame(score_rows)

                df.to_csv(score_csv_path, index=False)

    # ----------------------
    # Summarise & Plot
    # ----------------------
    results_csv = output_dir / "multi_pocket_docking_summary.csv"

    all_scores = []
    for pocket_dir in output_dir.glob("pocket_*"):
        score_csv = pocket_dir / "docking_scores.csv"
        if score_csv.exists() and score_csv.stat().st_size > 0:
            df = pd.read_csv(score_csv)
            all_scores.append(df)

    if all_scores:
        summary_df = pd.concat(all_scores, ignore_index=True)
        summary_df.to_csv(results_csv, index=False)
        logger.info(f"Saved docking summary to: {results_csv}")
    else:
        logger.warning("No docking scores found — summary CSV not created.")

    if config.get("plot_docking_results", True):
        plot_multi_pocket_scores(results_csv, output_dir)
    else:
        logger.info("Skipping result plotting.")

    logger.info("Multi-pocket docking workflow completed.")
