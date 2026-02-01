import yaml
import gc
from pathlib import Path

from workflows import register_workflow
from pipeline.logger import setup_logger
from pipeline.task_registry import get_task, list_tasks

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_workflow("train_model", description="Train and evaluate specified model")
def train_model(config_path: str):
    """
    Dynamically runs model training workflow as defined in a YAML config file.
    """
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    workflow_name = config.get("workflow_name", "train_pose_ranker")
    task_list = config.get("workflow", [])
    if not task_list:
        raise ValueError("No 'workflow' list found in config.")

    logger.info(f"[{workflow_name}] Starting workflow with tasks: {task_list}")

    current_data = {}

    for task_name in task_list:
        task_func = get_task(task_name)
        if task_func is None:
            available = ", ".join(sorted(list_tasks()))
            logger.error(f"❌ Task '{task_name}' not found in registry. Available tasks: {available}")
            raise ValueError(f"Task '{task_name}' not found in registry.")

        task_config = config.get(task_name, {})
        logger.info(f">>>>>>>>>> Running task: {task_name}")

        try:
            result = task_func(task_config, current_data)
        except Exception as e:
            logger.error(f"❌ Task '{task_name}' failed: {e}")
            raise

        if isinstance(result, dict):
            current_data.update(result)
        elif result is not None:
            current_data[task_name] = result

        gc.collect()

    logger.info(f"[{workflow_name}] Workflow complete.")
    return current_data

@register_workflow("vhh_active_learning", description="Run VHH generation active learning workflow")
def run_vhh_active_learning(config_path: str):
    """
    Dynamically runs the VHH active learning workflow as defined in a YAML config file.
    """
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    workflow_name = config.get("workflow_name", "vhh_active_learning")
    task_list = config.get("workflow", [])
    if not task_list:
        raise ValueError("No 'workflow' list found in config.")

    logger.info(f"[{workflow_name}] Starting workflow with tasks: {task_list}")

    current_data = {}

    for task_name in task_list:
        task_func = get_task(task_name)
        if task_func is None:
            available = ", ".join(sorted(list_tasks()))
            logger.error(f"❌ Task '{task_name}' not found in registry. Available tasks: {available}")
            raise ValueError(f"Task '{task_name}' not found in registry.")

        task_config = config.get(task_name, {})
        logger.info(f">>>>>>>>>> Running task: {task_name}")

        try:
            result = task_func(task_config, current_data)
        except Exception as e:
            logger.error(f"❌ Task '{task_name}' failed: {e}")
            raise

        if isinstance(result, dict):
            current_data.update(result)
        elif result is not None:
            current_data[task_name] = result

        gc.collect()

    logger.info(f"[{workflow_name}] Workflow complete.")
    return current_data