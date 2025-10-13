from workflows import register_workflow
from pipeline.logger import setup_logger
from pathlib import Path
from pipeline.task_registry import get_task

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_workflow(
    "pbpk_workflow",
    description="Predict PBPK parameters, assemble model, simulate, and analyse."
)
def run(config: dict):
    """
    Main workflow for PBPK simulation.
    Steps are defined in YAML: parameter prediction -> model assembly -> simulation -> analysis.
    """
    output_dir = Path(config.get("output_dir", "outputs/pbpk"))
    output_dir.mkdir(parents=True, exist_ok=True)
    config["output_dir"] = output_dir

    # --- Normalise predictor section (for safety) ---
    pbpk_section = config.setdefault("pbpk_parameter_prediction", {})
    predictors = pbpk_section.get("predictors", {})
    if not isinstance(predictors, dict):
        logger.warning("Predictors section malformed or missing — resetting to empty dict.")
        predictors = {}
    pbpk_section["predictors"] = predictors
    logger.info(f"[Config check] Predictors in YAML: {predictors}")

    # --- Workflow steps ---
    context = {}
    workflow_steps = config.get("workflow_serial", [])
    if not workflow_steps:
        raise ValueError("No workflow steps defined under 'workflow_serial' in config YAML.")

    for step in workflow_steps:
        logger.info(f"Running PBPK task: {step}")
        task_func = get_task(step)
        if not task_func:
            raise ValueError(f"Unknown task: {step}")

        # Log any per-step config sections for transparency
        step_section = config.get(step, None)
        if step_section:
            logger.debug(f"Subconfig for step '{step}': {step_section}")

        # Run the step
        result = task_func(config, **context)

            # --- Carry all returned data forward ---
        if isinstance(result, dict):
            # Keep all previous keys, update with any new ones from the task
            context.update({k: v for k, v in result.items() if v is not None})
        else:
            logger.debug(f"Task {step} returned non-dict result; context unchanged.")

    logger.info("PBPK workflow complete.")
    return context
