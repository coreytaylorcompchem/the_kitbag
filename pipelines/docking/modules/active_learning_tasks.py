import pandas as pd
import random
from pathlib import Path
from sklearn.ensemble import RandomForestRegressor
from lightgbm import LGBMRegressor
from rdkit import Chem
import numpy as np

from pipeline.task_registry import register_task, get_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)


@register_task("active_learn_docking", category="Active Learning", description="Active learning loop over docking.")
def active_learn_docking(backend, ligands, config, **kwargs):
    al_config = config['active_learning']
    n_initial = al_config['n_initial']
    batch_size = al_config['batch_size']
    n_cycles = al_config['n_cycles']
    seed = al_config.get('seed', 42)
    model_name = al_config.get('model', 'random_forest').lower()
    sample_size = al_config.get('sample_size', 5000)

    random.seed(seed)
    np.random.seed(seed)

    all_indices = list(range(len(ligands)))
    labeled_indices = random.sample(all_indices, n_initial)
    unlabeled_indices = list(set(all_indices) - set(labeled_indices))
    scores = {}

    use_confs = config.get('options', {}).get('use_conformer_generation', True)

    for cycle in range(n_cycles):
        logger.info(f">>>>>>>>>> Active Learning Cycle {cycle + 1}/{n_cycles} <<<<<<<<<<")

        # Select batch ligands to dock (those not already scored)
        current_batch = [ligands[i] for i in labeled_indices if ligands[i]['name'] not in scores]

        for ligand in current_batch:
            try:
                get_task("standardise_ligand")(backend, ligand, config)

                if use_confs:
                    get_task("generate_conformers")(backend, ligand, config)
                    get_task("cluster_conformers")(backend, ligand, config)

                get_task("convert_to_pdbqt")(backend, ligand, config)
                get_task("dock")(backend, ligand, config)

                scores[ligand['name']] = ligand.get('score')
                logger.info(f"Docked ligand {ligand['name']} with score: {ligand.get('score')}")

            except Exception as e:
                logger.error(f"Failed to prep/dock ligand {ligand['name']}: {e}")
                ligand['score'] = None
                scores[ligand['name']] = None

        # Train model on labeled data
        X, y = [], []
        for idx in labeled_indices:
            lig = ligands[idx]
            mol = Chem.MolFromSmiles(lig['smiles'])
            if mol is None or lig['name'] not in scores or scores[lig['name']] is None:
                continue
            fp = Chem.RDKFingerprint(mol)
            X.append(np.array(fp))
            y.append(scores[lig['name']])

        if not X:
            logger.warning("No valid labeled data to train model. Ending active learning.")
            break

        X = np.array(X)
        y = np.array(y)

        if model_name == "random_forest":
            model = RandomForestRegressor(random_state=seed, n_estimators=100)
            model.fit(X, y)

        elif model_name == "lightgbm":
            model = LGBMRegressor(random_state=seed, n_estimators=100)
            model.fit(X, y)

        else:
            raise ValueError(f"Unknown model specified in config: {model_name}")

        # --- Subsample unlabeled pool ---
        sample_size = min(sample_size, len(unlabeled_indices))
        sampled_unlabeled_indices = random.sample(unlabeled_indices, sample_size)

        X_unlabeled = []
        unlabeled_map = {}

        for idx in sampled_unlabeled_indices:
            mol = Chem.MolFromSmiles(ligands[idx]['smiles'])
            if mol is None:
                continue
            fp = Chem.RDKFingerprint(mol)
            X_unlabeled.append(np.array(fp))
            unlabeled_map[len(X_unlabeled) - 1] = idx

        if not X_unlabeled:
            logger.info("No valid unlabeled ligands remaining. Stopping.")
            break

        X_unlabeled = np.array(X_unlabeled)

        # Predict + estimate uncertainty
        if model_name == "random_forest":
            preds = model.predict(X_unlabeled)
            uncertainty = np.std([t.predict(X_unlabeled) for t in model.estimators_], axis=0)

        elif model_name == "lightgbm":
            preds_list = []
            for i in range(5):  # Ensemble of 5
                lgb_model = LGBMRegressor(random_state=seed + i, n_estimators=100)
                lgb_model.fit(X, y)
                preds_list.append(lgb_model.predict(X_unlabeled))
            preds_array = np.array(preds_list)
            preds = np.mean(preds_array, axis=0)
            uncertainty = np.std(preds_array, axis=0)

        # Acquisition: pick top uncertain
        acquisition = uncertainty
        selected_idx = np.argsort(-acquisition)[:batch_size]
        new_indices = [unlabeled_map[i] for i in selected_idx]

        labeled_indices.extend(new_indices)
        unlabeled_indices = list(set(unlabeled_indices) - set(new_indices))

        logger.info(f"Selected {len(new_indices)} new ligands for docking.")

        if not unlabeled_indices:
            logger.info("No more unlabeled ligands remaining. Ending active learning.")
            break

    # Save scores
    df = pd.DataFrame({'name': list(scores.keys()), 'score': list(scores.values())})
    df.to_csv(config['output_dir'] / "docking_scores.csv", index=False)
