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

    # Use a single context for all tasks
    current_data = {}
    # Inject full YAML for logging
    current_data["full_config"] = config

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

@register_workflow("run_adme_inference", description="Run ADME model inference on SMILES CSV")
def run_adme_inference(config_path: str):
    """
    Dynamically runs ADME inference workflow as defined in a YAML config file.
    """

    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    workflow_name = config.get("workflow_name", "inference_adme_model")
    task_list = config.get("workflow", [])

    if not task_list:
        raise ValueError("No 'workflow' list found in config.")

    logger.info(f"[{workflow_name}] Starting inference workflow with tasks: {task_list}")

    current_data = {}

    # Optional but useful for debugging / reproducibility
    current_data["full_config"] = config

    for task_name in task_list:
        task_func = get_task(task_name)

        if task_func is None:
            available = ", ".join(sorted(list_tasks()))
            logger.error(
                f"❌ Task '{task_name}' not found in registry. Available tasks: {available}"
            )
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

    logger.info(f"[{workflow_name}] Inference workflow complete.")
    return current_data

@register_workflow(
    "migrate_adme_models_to_mlflow",
    description="Migrate existing ADME models and evaluation artifacts to MLflow",
)
def migrate_adme_models_to_mlflow(config_path: str):
    """
    Run the historical ADME model migration workflow defined in a YAML file.

    This workflow uploads existing checkpoints, preprocessing artifacts,
    training configurations, evaluation metrics, and plots to MLflow without
    retraining the models.
    """
    config_path = Path(config_path)

    if not config_path.exists():
        raise FileNotFoundError(
            f"Config file not found: {config_path}"
        )

    with open(config_path, "r", encoding="utf-8") as file:
        config = yaml.safe_load(file)

    if not isinstance(config, dict):
        raise ValueError(
            f"Workflow configuration must be a YAML mapping: {config_path}"
        )

    workflow_name = config.get(
        "workflow_name",
        "migrate_adme_models_to_mlflow",
    )

    task_list = config.get("workflow", [])

    if not task_list:
        raise ValueError(
            "No 'workflow' list found in config."
        )

    logger.info(
        f"[{workflow_name}] Starting MLflow migration workflow "
        f"with tasks: {task_list}"
    )

    current_data = {
        "full_config": config,
        "config_path": str(config_path.resolve()),
    }

    for task_name in task_list:
        task_func = get_task(task_name)

        if task_func is None:
            available = ", ".join(
                sorted(list_tasks())
            )

            logger.error(
                f"❌ Task '{task_name}' not found in registry. "
                f"Available tasks: {available}"
            )

            raise ValueError(
                f"Task '{task_name}' not found in registry."
            )

        task_config = config.get(task_name, {})

        logger.info(
            f">>>>>>>>>> Running task: {task_name}"
        )

        try:
            result = task_func(
                task_config,
                current_data,
            )

        except Exception as error:
            logger.error(
                f"❌ Task '{task_name}' failed: {error}"
            )
            raise

        if isinstance(result, dict):
            current_data.update(result)

        elif result is not None:
            current_data[task_name] = result

        gc.collect()

    logger.info(
        f"[{workflow_name}] MLflow migration workflow complete."
    )

    return current_data