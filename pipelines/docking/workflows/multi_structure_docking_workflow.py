import yaml
import csv
import pandas as pd
from pathlib import Path
from rdkit import Chem

from backends import get_backend
from modules.protein_preparation import ProteinPreparer
from pipeline.task_registry import get_task
from workflows import register_workflow

from workflows.utils.docking import (
    generate_ligands_csv_from_txt,
    summarise_docking_results,
    validate_config,
    plot_multi_structure_scores,
    extract_crystal_ligand_center,
    extract_scores_from_docking_output
)

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_workflow("multi_structure_docking", description="Dock ligands into multiple PDB structures.")
def run(config_path: str):
    import shutil

    # ----------------------
    # Load YAML Config
    # ----------------------
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

    docking_cfg = config["docking"]
    output_dir = Path(docking_cfg["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    config["output_dir"] = output_dir

    # ----------------------
    # Ligands
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

    ligands = []
    with open(ligands_csv_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ligands.append({'name': row['name'], 'smiles': row['smiles']})

    # ----------------------
    # Backend
    # ----------------------
    backend_name = config['backend']['name']
    backend_kwargs = {k: v for k, v in config['backend'].items() if k != 'name'}
    backend = get_backend(backend_name, **backend_kwargs)

    # ----------------------
    # Multi-structure Docking
    # ----------------------
    protein_dir = Path(config.get('protein', {}).get('protein_directory', None))
    if not protein_dir.exists():
        raise FileNotFoundError(f"Protein directory not found: {protein_dir}")

    pdb_files = sorted(protein_dir.glob("*.pdb"))
    if not pdb_files:
        raise FileNotFoundError(f"No PDB files found in: {protein_dir}")

    logger.info(f"Found {len(pdb_files)} PDB files for docking.")

    for i, pdb_path in enumerate(pdb_files):
        pdb_path = Path(pdb_path)
        logger.info(f"\n>>>>>>>>>> Docking to structure {i+1}/{len(pdb_files)}: {pdb_path.name} <<<<<<<<<<")

        # Setup per-protein subdir
        protein_output_dir = output_dir / pdb_path.stem
        protein_output_dir.mkdir(parents=True, exist_ok=True)
        config["output_dir"] = protein_output_dir  # update config for this round

        # Prepare protein
        protein_preparer = ProteinPreparer(
            pdb_path=pdb_path,
            work_dir=protein_output_dir,
            pH=config.get('protein', {}).get('pH', 7.4)
        )
        receptor_pdbqt = protein_preparer.prepare()
        backend.cache["receptor_pdbqt"] = receptor_pdbqt

        # Determine binding site center from crystal ligand
        center = extract_crystal_ligand_center(pdb_path)
        config['docking']['center'] = center

        # Ensure config['protein']['pdb_path'] is set
        config.setdefault("protein", {})["pdb_path"] = str(pdb_path)

        for ligand in ligands:
            logger.info(f"Processing ligand: {ligand['name']} against {pdb_path.name}")
            docking_outputs = []

            for step in config.get("workflow", []):
                task_func = get_task(step)
                if not task_func:
                    raise ValueError(f"Workflow step '{step}' is not a registered task.")

                result = task_func(
                    backend,
                    ligand,
                    config,
                    protein_id=pdb_path.stem,
                    pdb_path=pdb_path
                )

                if isinstance(result, list):
                    docking_outputs.extend(result)
                elif result:
                    docking_outputs.append(result)

            logger.debug(f"docking_outputs = {docking_outputs}")

            if docking_outputs:
                pocket_dir = protein_output_dir / "pocket_0"
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
                                "structure": pdb_path.stem
                            })
                
                # Check there is valid docking output
                logger.debug(score_rows)
                
                score_csv_path = pocket_dir / "docking_scores.csv"
                
                # Check if csv is already present 

                if score_csv_path.exists():
                    if score_csv_path.stat().st_size > 0:
                        existing_df = pd.read_csv(score_csv_path)
                        df = pd.concat([existing_df, pd.DataFrame(score_rows)], ignore_index=True)
                    else:
                        logger.warning(f"Existing docking_scores.csv is empty — will overwrite it.")
                        df = pd.DataFrame(score_rows)
                else:
                    df = pd.DataFrame(score_rows)

                df.to_csv(score_csv_path, index=False)



    # ----------------------
    # Summarise
    # ----------------------
    final_summary_csv = output_dir / "multi_structure_docking_summary.csv"
    summarise_docking_results(output_dir, final_summary_csv)

    if config.get("plot_summarised_docking_results", True):
        plot_multi_structure_scores(final_summary_csv, output_dir)
    else:
        logger.info("Skipping result plotting.")

    logger.info("Multi-structure docking workflow completed.")
