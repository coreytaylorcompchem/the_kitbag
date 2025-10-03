import yaml
import csv
from pathlib import Path
from workflows import register_workflow
from modules.protein_preparation import ProteinPreparer
from pipeline.task_registry import get_task
from backends import get_backend
from workflows.utils.docking import (
    generate_ligands_csv_from_txt,
    get_docking_box,
    extract_scores_from_docking_output,
)
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_workflow("active_learning_docking", description="prepare, dock and score with Active Learning loop.")
def run(config_path: str):
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    output_dir = Path(config['docking']['output_dir'])
    output_dir.mkdir(parents=True, exist_ok=True)
    config['output_dir'] = output_dir

    # Ligand loading
    ligands_txt_path = Path(config.get('ligands_txt', 'ligands.txt'))
    ligands_csv_path = output_dir / 'ligands.csv'
    generate_ligands_csv_from_txt(ligands_txt_path, ligands_csv_path)

    ligands = []
    with open(ligands_csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            ligands.append({'name': row['name'], 'smiles': row['smiles']})

    logger.info(f"Loaded {len(ligands)} ligands.")

    # Protein + backend
    protein_preparer = ProteinPreparer(config['protein']['pdb_path'], output_dir, pH=config['protein']['pH'])
    receptor_pdbqt = protein_preparer.prepare()
    config['receptor_pdbqt'] = receptor_pdbqt

    backend = get_backend(config['backend']['name'], **{k: v for k, v in config['backend'].items() if k != 'name'})
    backend.cache['receptor_pdbqt'] = receptor_pdbqt

    # Determine docking box center + size if not already set
    if 'center' not in config['docking'] or 'size' not in config['docking']:
        center, size = get_docking_box(config, protein_preparer)
        config['docking']['center'] = center
        config['docking']['size'] = size
        logger.info(f"Auto-detected docking box: center={center}, size={size}")

    # Enforce conformer generation config flag
    if 'options' not in config:
        config['options'] = {}
    config['options']['use_conformer_generation'] = config['options'].get('use_conformer_generation', True)

    # --- PATCH: wrap backend.dock to handle score extraction ---
    original_dock = backend.dock

    def patched_dock(ligand, config, **kwargs):
        result = original_dock(ligand, config, **kwargs)

        # Ensure result is a Path (handle list/dict/str)
        if isinstance(result, list):
            result_path = Path(result[0])
        elif isinstance(result, (str, Path)):
            result_path = Path(result)
        else:
            raise TypeError(f"Unexpected docking result type: {type(result)}")

        ligand['docked_output'] = str(result_path)

        try:
            scores = extract_scores_from_docking_output(result_path)
            if scores:
                ligand['score'] = min(s['score'] for s in scores)
                ligand['poses'] = scores
            else:
                logger.warning(f"No docking scores extracted for ligand: {ligand['name']}")
                ligand['score'] = None
        except Exception as e:
            logger.warning(f"Failed to extract scores for ligand {ligand['name']}: {e}")
            ligand['score'] = None

        return result_path

    backend.dock = patched_dock
    # ----------------------------------------------------------

    # Call AL task
    al_task = get_task("active_learn_docking")
    if al_task is None:
        raise RuntimeError("Active learning docking task 'active_learn_docking' not found.")

    al_task(backend, ligands, config)

    # Restore original dock method
    backend.dock = original_dock

    # Scores are saved by the AL task, but optionally save summary here again if needed
    scores_csv = output_dir / "docking_scores.csv"
    if scores_csv.exists():
        logger.info(f"Active learning docking scores saved to: {scores_csv}")
    else:
        logger.warning(f"No docking scores file found at {scores_csv}")
    
    logger.info("Active learning docking completed.")
