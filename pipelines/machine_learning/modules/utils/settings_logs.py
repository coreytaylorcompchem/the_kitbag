from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def log_pipeline_config(config, context):
    """
    Logs key configuration and context values at the start of the VHH active learning workflow.
    """

    workflow_name = config.get("workflow_name", "VHH_active_learning")
    workflow_steps = config.get("workflow", [])

    logger.info("-" * 60)
    logger.info(f"Workflow: {workflow_name}")
    logger.info(f"Workflow steps: {workflow_steps}")

    # ----------------------
    # Dataset settings
    # ----------------------
    ds_cfg = config.get("load_seed_dataset", {})
    logger.info(f"CSV path: {ds_cfg.get('csv_path', 'Unknown')}")
    logger.info(f"Expected sequence length: {ds_cfg.get('expected_len', 'Unknown')}")
    logger.info(f"Multi-condition properties: {ds_cfg.get('multi_condition_props', [])}")

    # ----------------------
    # ESM model settings
    # ----------------------
    esm_cfg = config.get("load_esm_model", {})
    logger.info(f"ESM device: {esm_cfg.get('device', 'cpu')}")
    logger.info(f"Pretrained model: {esm_cfg.get('pretrained_model', 'esm1b_t33_650M_UR50S')}")
    logger.info(f"Batch size: {esm_cfg.get('batch_size', 32)}")

    # ----------------------
    # Active learning settings
    # ----------------------
    al_cfg = config.get("active_learning_rounds", {})
    logger.info(f"Number of rounds: {al_cfg.get('n_rounds', 3)}")
    logger.info(f"Samples per seed: {al_cfg.get('samples_per_seed', 100)}")
    logger.info(f"Top candidates per round: {al_cfg.get('top_candidates_per_round', 1000)}")
    logger.info(f"Model type: CatBoostRegressor")
    logger.info(f"Data directory: {al_cfg.get('data_dir', 'data_vhh')}")
    logger.info(f"Plots directory: {al_cfg.get('plots_dir', 'plots_vhh')}")

    # Optional debug info
    logger.debug(f"CatBoost params: {al_cfg.get('catboost_params', {})}")
    logger.debug(f"PCA components: {al_cfg.get('n_components', 512)}")
    logger.debug(f"UMAP settings: n_neighbors={al_cfg.get('umap_n_neighbors', 15)}, "
                 f"min_dist={al_cfg.get('umap_min_dist', 0.1)}, n_epochs={al_cfg.get('umap_n_epochs', 50)}")

    # Multi-condition properties from context
    logger.debug(f"Multi-condition properties in context: {context.get('multi_condition_props', [])}")

    logger.info("-" * 60)
