from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def log_pipeline_config(config, context, full_config=None):
    logger.info("-" * 60)

    # Workflow-level info
    if full_config:
        logger.info(f"Workflow: {full_config.get('workflow_name', 'VHH_active_learning')}")
        logger.info(f"Workflow steps: {full_config.get('workflow', [])}")

    # Task-level info
    keys_to_log = [
        "csv_path",
        "expected_len",
        "multi_condition_props",
        "n_rounds",
        "samples_per_seed",
        "top_candidates_per_round",
        "min_per_ancestor",
        "noise_scale",
        "batch_size",
        "n_model_ensemble",
        "max_mutants",
        "data_dir",
        "plots_dir",
    ]
    for key in keys_to_log:
        value = config.get(key)
        if value is None and full_config:
            # fallback to task slice in full_config
            value = full_config.get("active_learning_rounds", {}).get(key)
        logger.info(f"{key.replace('_', ' ').capitalize()}: {value}")

    # Model info
    catboost_params = config.get("catboost_params")
    if catboost_params is None and full_config:
        catboost_params = full_config.get("active_learning_rounds", {}).get("catboost_params")
    logger.info(f"Model type: CatBoostRegressor")
    logger.debug(f"CatBoost params: {catboost_params}")

    # Device info
    logger.info(f"ESM device: {context.get('device', 'Unknown')}")
    logger.info(f"Pretrained model: {context.get('esm_model', 'Unknown')}")

    logger.info("-" * 60)
