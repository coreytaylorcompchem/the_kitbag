from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def log_pipeline_config(config, context, full_config=None):
    """
    Logs key configuration and context values at the start of the VHH active learning workflow.
    
    Args:
        config (dict): Task-specific YAML slice.
        context (dict): Runtime context.
        full_config (dict, optional): Full YAML dict, for workflow-level logging.
    """
    logger.info("-" * 60)

    if full_config:
        logger.info(f"Workflow: {full_config.get('workflow_name', 'VHH_active_learning')}")
        logger.info(f"Workflow steps: {full_config.get('workflow', [])}")

    # Task-specific info
    logger.info(f"CSV path: {config.get('csv_path', 'Unknown')}")
    logger.info(f"Expected sequence length: {config.get('expected_len', 'Unknown')}")
    logger.info(f"Multi-condition properties: {config.get('multi_condition_props', [])}")

    # Model / training info
    logger.info(f"Model type: CatBoostRegressor")
    logger.info(f"Data directory: {config.get('data_dir', 'data_vhh')}")
    logger.info(f"Plots directory: {config.get('plots_dir', 'plots_vhh')}")

    logger.info("-" * 60)
