from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def log_pipeline_config(config, context, full_config=None):
    """
    Logs key configuration and context values at the start of the VHH active learning workflow.

    Args:
        config (dict): Task-specific YAML slice (load_seed_dataset).
        context (dict): Runtime context.
        full_config (dict, optional): Entire YAML for workflow-level logging.
    """
    workflow_name = full_config.get("workflow_name", "VHH_active_learning") if full_config else config.get("workflow_name", "VHH_active_learning")
    workflow_steps = full_config.get("workflow", []) if full_config else []

    logger.info("-" * 60)
    logger.info(f"Workflow: {workflow_name}")
    logger.info(f"Workflow steps: {workflow_steps}")

    # Task-specific info
    logger.info(f"CSV path: {config.get('csv_path', 'Unknown')}")
    logger.info(f"Expected len: {config.get('expected_len', 'Unknown')}")
    logger.info(f"Multi condition props: {config.get('multi_condition_props', [])}")

    # Now workflow-wide / active_learning_rounds keys
    al_config = full_config.get("active_learning_rounds", {}) if full_config else {}
    logger.info(f"N rounds: {al_config.get('n_rounds', 'Unknown')}")
    logger.info(f"Samples per seed: {al_config.get('samples_per_seed', 'Unknown')}")
    logger.info(f"Top candidates per round: {al_config.get('top_candidates_per_round', 'Unknown')}")
    logger.info(f"Min per ancestor: {al_config.get('min_per_ancestor', 'Unknown')}")
    logger.info(f"Noise scale: {al_config.get('noise_scale', 'Unknown')}")
    logger.info(f"Batch size: {al_config.get('batch_size', 'Unknown')}")
    logger.info(f"N model ensemble: {al_config.get('n_model_ensemble', 'Unknown')}")
    logger.info(f"Max mutants: {al_config.get('max_mutants', 'Unknown')}")
    logger.info(f"Data dir: {al_config.get('data_dir', 'Unknown')}")
    logger.info(f"Plots dir: {al_config.get('plots_dir', 'Unknown')}")

    # Model / training info
    logger.info(f"Model type: CatBoostRegressor")
    esm_config = full_config.get("load_esm_model", {}) if full_config else {}
    logger.info(f"ESM device: {esm_config.get('device', 'Unknown')}")
    logger.info(f"Pretrained model: {esm_config.get('pretrained_model', 'Unknown')}")

    logger.info("-" * 60)
