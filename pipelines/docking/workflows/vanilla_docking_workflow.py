import os
import yaml
import csv

import pandas as pd
from pathlib import Path

from concurrent.futures import ThreadPoolExecutor, as_completed

from backends import get_backend, discover_backends
from modules.protein_preparation import ProteinPreparer
from pipeline.task_registry import get_task
from workflows import register_workflow

from workflows.utils.docking import (
    generate_ligands_csv_from_txt,
    validate_config,
    extract_scores_from_docking_output,
    get_docking_box,
)

from workflows.utils.rdkit import compute_rdkit_descriptors

from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

def run_single_ligand(ligand, backend, config, workflow_steps):
    docking_outputs = []

    for step in workflow_steps:
        task_func = get_task(step)
        if not task_func:
            raise ValueError(f"Workflow step '{step}' is not a registered task.")

        result = task_func(backend, ligand, config)

        if step == "dock" and result:
            if isinstance(result, list):
                docking_outputs.extend(result)
            else:
                docking_outputs.append(result)

    # ----------------------
    # Extract docking scores
    # ----------------------
    all_scores = []

    for entry in docking_outputs:
        docked_outputs = entry.get("docked_output")

        if not docked_outputs:
            continue

        if not isinstance(docked_outputs, (list, tuple)):
            docked_outputs = [docked_outputs]

        for sdf_path in docked_outputs:
            scores = extract_scores_from_docking_output(sdf_path)

            for s in scores:
                s["sdf_path"] = str(sdf_path)

            all_scores.extend(scores)

    if all_scores:
        ligand["docking_scores"] = all_scores
        best = min(all_scores, key=lambda x: x["score"])
        ligand["docking_score"] = best["score"]
        ligand["best_pose_idx"] = best.get("pose_idx")
        ligand["best_pose_rank"] = best.get("pose_rank")

    return ligand["name"]


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
    backend_cfg = config.get("backend", {})

    logger.info("---------------------------------------------")
    logger.info("Loaded YAML configuration:")
    logger.info(f"Backend / using GPU: {backend_cfg.get('name')} / {backend_cfg.get('use_gpu')}")
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
        logger.info(f"Found ligand file - generating ligands.csv.")
        generate_ligands_csv_from_txt(ligands_txt_path, ligands_csv_path)
    elif ligands_csv_path.exists():
        logger.info(f"Found ligands.csv - using it directly.")
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
            lig = {
                "name": row["name"],
                "smiles": row["smiles"],
            }

            # RDKit descriptors (once per ligand)
            lig.update(compute_rdkit_descriptors(row["smiles"]))

            ligands.append(lig)

    # ----------------------
    # Execute Workflow Steps
    # ----------------------
    workflow_steps = config.get("workflow", [
        "standardise_ligand",
        "generate_conformers",
        "cluster_conformers",
        "save_final_conformers",
        "convert_to_pdbqt",
        "dock",
    ])

    # --------------------------------------
    # Separate per-ligand and dataset tasks
    # --------------------------------------
    DATASET_LEVEL_TASKS = {"post_processing"}

    per_ligand_steps = [
        step for step in workflow_steps
        if step not in DATASET_LEVEL_TASKS
    ]

    dataset_steps = [
        step for step in workflow_steps
        if step in DATASET_LEVEL_TASKS
    ]

    # Handle optional conformer generation
    options = config.get("options", {})
    use_conformer_generation = options.get("use_conformer_generation", True)

    if not use_conformer_generation:
        logger.warning("Conformer generation steps disabled - skipping.")
        workflow_steps = [
            step for step in workflow_steps
            if step not in ("generate_conformers", "cluster_conformers", "save_final_conformers")
        ]

    # ----------------------
    # Parallel code (per ligand)
    # ----------------------
    n_cores = os.cpu_count() or 20
    n_workers = max(1, n_cores - 2)  # reserve 2 cores
    logger.info(f"Using {n_workers} parallel workers (out of {n_cores} cores).")

    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = [
            executor.submit(run_single_ligand, lig, backend, config, per_ligand_steps)
            for lig in ligands
        ]
        for future in as_completed(futures):
            try:
                ligand_name = future.result()
            except Exception as e:
                logger.error(f"❌ Error during ligand processing: {e}")
        
        results = []

        LIGAND_EXCLUDE_KEYS = {
            "docking_scores",
            "best_pose_idx",
            "best_pose_rank",
            "adme",
            "pdbqt_paths",
            "pdbqt_path"
        }

        for lig in ligands:
            base = {
                k: v for k, v in lig.items()
                if k not in LIGAND_EXCLUDE_KEYS
            }

            # ADME (copied to every pose)
            if "adme" in lig:
                base.update(lig["adme"])

            docking_scores = lig.get("docking_scores", [])

            if docking_scores:
                for entry in docking_scores:
                    row = base.copy()
                    row.update({
                        "pose_idx": entry.get("pose_idx"),
                        "pose_rank": entry.get("pose_rank"),
                        "docking_score": entry.get("score"),
                        "sdf_path": entry.get("sdf_path"),
                    })
                    results.append(row)
            else:
                # ligand with no docking result
                results.append(base)

        df = pd.DataFrame(results)

        out_csv = output_dir / "docking_results.csv"
        df.to_csv(out_csv, index=False)

        logger.info(f"Saved combined ADME/docking results to: {out_csv}")

    # --------------------------------------
    # Run dataset-level workflow steps
    # --------------------------------------
    for step in dataset_steps:
        task_func = get_task(step)
        if not task_func:
            raise ValueError(
                f"Workflow step '{step}' is not a registered task."
            )

        logger.info(f"Running dataset-level task: {step}")

        try:
            task_func(backend, None, config)
        except Exception as e:
            logger.error(f"❌ Dataset task '{step}' failed: {e}")

    logger.info("Vanilla docking workflow completed.")
