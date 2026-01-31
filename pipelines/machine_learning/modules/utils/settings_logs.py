from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def log_pipeline_config(config, context):
    """
    Logs key configuration and context values at the start of the VHH active learning workflow.
    
    Args:
        config (dict): YAML configuration for the workflow.
        context (dict): Runtime context including multi-condition props, device, etc.
    """
    workflow_name = config.get("workflow_name", "VHH_active_learning")
    workflow_steps = config.get("workflow", [])
    
    logger.info("-" * 60)
    
    # General workflow info
    logger.info(f"Number of rounds: {config.get('n_rounds', 3)}")
    logger.info(f"Samples per seed: {config.get('samples_per_seed', 100)}")
    logger.debug(f"Max mutants per sequence: {config.get('max_mutants', 3)}")
    logger.info(f"Top candidates per round: {config.get('top_candidates_per_round', 1000)}")
    logger.debug(f"Batch size for embeddings / predictions: {config.get('batch_size', 32)}")
    logger.debug(f"Random seed (if any): {config.get('random_seed', 'None')}")
    
    # Multi-condition and sequence info
    logger.info(f"Multi-condition properties: {context.get('multi_condition_props', [])}")
    logger.debug(f"Expected sequence length: {config.get('expected_len', 'Unknown')}")
    mutable_positions = context.get("mutable_positions")
    if mutable_positions is not None:
        logger.debug(f"Number of mutable positions: {len(mutable_positions)}")
    logger.debug(f"Number of initial seeds: {len(context.get('current_seeds', []))}")
    
    # Model / training info
    logger.info(f"Model type: CatBoostRegressor")
    logger.debug(f"CatBoost params: {config.get('catboost_params', {})}")
    device = context.get("device")
    logger.debug(f"Device for ESM embeddings / model: {device}")
    
    # Visualization info
    logger.info(f"Data directory: {config.get('data_dir', 'data_vhh')}")
    logger.info(f"Plots directory: {config.get('plots_dir', 'plots_vhh')}")
    logger.debug(f"UMAP settings: n_neighbors={config.get('umap_n_neighbors', 15)}, "
                f"min_dist={config.get('umap_min_dist', 0.1)}, n_epochs={config.get('umap_n_epochs', 50)}")
    
    logger.info("-" * 60)
