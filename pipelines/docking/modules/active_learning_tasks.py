import pandas as pd
import random
from pathlib import Path
from sklearn.ensemble import RandomForestRegressor
import lightgbm as lgb
from rdkit import Chem
import numpy as np

from modules.utils.plot_performance import plot_performance
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
    sample_size = al_config.get('sample_size', None)  # Optional sampling
    model_name = al_config.get('model', 'random_forest')
    acquisition_method = al_config.get('acquisition', 'uncertainty_sampling')
    run_baseline = al_config.get('run_random_baseline', False)

    random.seed(seed)
    np.random.seed(seed)

    all_indices = list(range(len(ligands)))
    labeled_indices = random.sample(all_indices, n_initial)
    unlabeled_indices = list(set(all_indices) - set(labeled_indices))
    scores = {}

    use_confs = config.get('options', {}).get('use_conformer_generation', True)

    # DataFrames to track per-cycle results
    al_results = []
    random_baseline_results = []

    def fingerprint_array(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        fp = Chem.RDKFingerprint(mol)
        return np.array(fp)

    def select_acquisition(X_unlabeled, model, method):
        if method == "uncertainty_sampling":
            # Uncertainty = std deviation of estimators predictions (RF or LGB model with n_estimators)
            if hasattr(model, 'estimators_'):  # RF
                preds_per_tree = np.array([est.predict(X_unlabeled) for est in model.estimators_])
                uncertainty = np.std(preds_per_tree, axis=0)
            elif hasattr(model, 'booster_'):  # LightGBM model
                # LightGBM does not expose individual tree predictions easily; fallback to predictive variance approx
                preds = model.predict(X_unlabeled)
                uncertainty = np.abs(preds - np.mean(preds))  # crude proxy, improve if desired
            else:
                uncertainty = np.zeros(len(X_unlabeled))
            return uncertainty

        elif method == "greedy_max_score":
            preds = model.predict(X_unlabeled)
            # If lower score is better, pick lowest predicted score (highest negative)
            return -preds

        elif method == "random_sampling":
            return np.random.rand(len(X_unlabeled))

        else:
            logger.warning(f"Unknown acquisition method {method}, defaulting to random.")
            return np.random.rand(len(X_unlabeled))

    for cycle in range(n_cycles):
        logger.info(f">>>>>>>>>> Active Learning Cycle {cycle + 1}/{n_cycles} <<<<<<<<<<")

        # Dock all ligands in labeled_indices that have not been scored yet
        current_batch = [ligands[i] for i in labeled_indices if ligands[i]['name'] not in scores]

        for ligand in current_batch:
            try:
                get_task("standardise_ligand")(backend, ligand, config)

                if use_confs:
                    get_task("generate_conformers")(backend, ligand, config)
                    get_task("cluster_conformers")(backend, ligand, config)

                get_task("convert_to_pdbqt")(backend, ligand, config)
                get_task("dock")(backend, ligand, config)

                ligand_score = ligand.get('score')
                scores[ligand['name']] = ligand_score
                logger.info(f"Docked ligand {ligand['name']} with score: {ligand_score}")

                # Track results for plotting
                al_results.append({
                    'cycle': cycle + 1,
                    'name': ligand['name'],
                    'smiles': ligand['smiles'],
                    'score': ligand_score
                })
            except Exception as e:
                logger.error(f"Failed to prep/dock ligand {ligand['name']}: {e}")
                ligand['score'] = None
                scores[ligand['name']] = None

        # Prepare training data (only ligands with valid scores)
        X, y = [], []
        for idx in labeled_indices:
            lig = ligands[idx]
            mol_fp = fingerprint_array(lig['smiles'])
            if mol_fp is None or lig['name'] not in scores or scores[lig['name']] is None:
                continue
            X.append(mol_fp)
            y.append(scores[lig['name']])

        if not X:
            logger.warning("No valid labeled data to train model. Ending active learning.")
            break

        X = np.array(X)
        y = np.array(y)

        # Train model
        if model_name == "random_forest":
            model = RandomForestRegressor(random_state=seed)
            model.fit(X, y)
        elif model_name == "lightgbm":
            train_data = lgb.Dataset(X, label=y)
            params = {'objective': 'regression', 'verbose': -1, 'seed': seed}
            model = lgb.train(params, train_data, num_boost_round=100)
        else:
            raise ValueError(f"Unsupported model {model_name}")

        # Subsample unlabeled set if sample_size is set
        if sample_size is not None and len(unlabeled_indices) > sample_size:
            sampled_unlabeled = random.sample(unlabeled_indices, sample_size)
        else:
            sampled_unlabeled = unlabeled_indices

        # Prepare fingerprints for unlabeled set
        X_unlabeled = []
        unlabeled_map = {}
        for idx_pos, idx in enumerate(sampled_unlabeled):
            mol_fp = fingerprint_array(ligands[idx]['smiles'])
            if mol_fp is None:
                continue
            X_unlabeled.append(mol_fp)
            unlabeled_map[len(X_unlabeled) - 1] = idx

        if not X_unlabeled:
            logger.info("No valid unlabeled ligands remaining. Stopping.")
            break

        X_unlabeled = np.array(X_unlabeled)

        # Acquisition
        acquisition_scores = select_acquisition(X_unlabeled, model, acquisition_method)

        # Select top batch_size indices by acquisition score
        selected_idx = np.argsort(-acquisition_scores)[:batch_size]
        new_indices = [unlabeled_map[i] for i in selected_idx]

        # Add selected ligands to labeled set and remove from unlabeled
        labeled_indices.extend(new_indices)
        unlabeled_indices = list(set(unlabeled_indices) - set(new_indices))

        logger.info(f"Selected {len(new_indices)} new ligands for docking.")

        # --- Random baseline (randomly select same batch size from unlabeled) ---
        if run_baseline and len(unlabeled_indices) >= batch_size:
            random_sample_indices = random.sample(unlabeled_indices, batch_size)
            for idx in random_sample_indices:
                ligand = ligands[idx]
                try:
                    get_task("standardise_ligand")(backend, ligand, config)

                    if use_confs:
                        get_task("generate_conformers")(backend, ligand, config)
                        get_task("cluster_conformers")(backend, ligand, config)

                    get_task("convert_to_pdbqt")(backend, ligand, config)
                    get_task("dock")(backend, ligand, config)

                    baseline_score = ligand.get('score')
                    if baseline_score is not None:
                        random_baseline_results.append({
                            'cycle': cycle + 1,
                            'name': ligand['name'],
                            'smiles': ligand['smiles'],
                            'score': baseline_score
                        })
                    logger.info(f"Random baseline docked ligand {ligand['name']} with score: {baseline_score}")
                except Exception as e:
                    logger.warning(f"Random baseline ligand {ligand['name']} failed: {e}")

        if len(unlabeled_indices) == 0:
            logger.info("No more unlabeled ligands remaining. Ending active learning.")
            break

    # Save AL scores per cycle
    cycle_df = pd.DataFrame(al_results)
    cycle_df.to_csv(config['output_dir'] / "per_cycle_scores.csv", index=False)

    # Save baseline scores if applicable
    if run_baseline and random_baseline_results:
        baseline_df = pd.DataFrame(random_baseline_results)
        baseline_df.to_csv(config['output_dir'] / "baseline_scores.csv", index=False)
    else:
        baseline_df = None

    # Plot results with batch and sample size in title
    plot_performance(
        per_cycle_df=cycle_df,
        output_dir=config['output_dir'],
        sample_size=sample_size,
        batch_size=batch_size,
        baseline_df=baseline_df
    )

    # Save overall docking scores (final)
    final_scores_df = pd.DataFrame({'name': list(scores.keys()), 'score': list(scores.values())})
    final_scores_df.to_csv(config['output_dir'] / "docking_scores.csv", index=False)
