from workflows import register_workflow, get_workflow
from pipeline.task_registry import get_task
from pipeline.backend_registry import get_backend_class

from pipeline.logger import setup_logger
logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_workflow(
    'staged_workflow',
    description="Performs single or multi-stage calculations."
)
def run_staged_workflow(xyz_file, config):
    steps = config.get('workflow', [])
    current_xyz = xyz_file
    last_result = {}

    for step in steps:
        logger.info(f">>> Running workflow step: {step} <<<")
        task_func = get_task(step)

        if task_func:
            step_config = config.get(step)
            if not step_config:
                raise ValueError(f"Missing config block for step: '{step}'")

            backend_name = step_config.get('backend')
            if not backend_name:
                raise ValueError(f"Missing 'backend' for step: {step}")
            
            BackendClass = get_backend_class(backend_name)
            backend = BackendClass(config)  # pass global config

            # For steps that require wavefunction input (like qtaim, mesp),
            # pass wfn_file path if available, else wfn object, else fallback to current_xyz
            if last_result and step in ("qtaim", "mesp"):
                if last_result.get("wfn_file"):
                    input_for_task = last_result["wfn_file"]
                elif last_result.get("wfn"):
                    input_for_task = last_result["wfn"]
                else:
                    logger.warning(f"No wavefunction found in last_result for step '{step}', passing current_xyz")
                    input_for_task = current_xyz

                # Copy method/basis from last_result if missing in step_config
                for key in ("method", "basis"):
                    if key not in step_config and key in last_result:
                        step_config[key] = last_result[key]
            else:
                input_for_task = current_xyz

            logger.debug(f"Passing to task '{step}': {type(input_for_task)} - {input_for_task}")

            result = task_func(backend, input_for_task, step_config, config)
            logger.debug(f"Result from step '{step}': {result}")

            # Update last_result and current_xyz depending on step and result
            if step in ('single_point', 'optimise'):
                if isinstance(result, dict):
                    last_result = result

                    # Prefer to update current_xyz with wfn_file if available
                    if "wfn_file" in result and result["wfn_file"]:
                        current_xyz = result["wfn_file"]
                        logger.debug(f"Updated current_xyz to wfn_file: {current_xyz}")
                    elif "wfn" in result and result["wfn"]:
                        current_xyz = result["wfn"].molecule().save_string_xyz()
                        logger.debug("Updated current_xyz to xyz string from wavefunction object.")
                    elif "geometry" in result:
                        current_xyz = result["geometry"]
                        logger.debug("Updated current_xyz to geometry string.")
                    else:
                        logger.debug("No update to current_xyz performed; using previous value.")
            elif getattr(task_func, 'modifies_geometry', False) and isinstance(result, str):
                current_xyz = result
                last_result = {"geometry": result}
                logger.debug(f"Task '{step}' modified geometry; updated current_xyz.")
            else:
                # For other steps, update last_result if result is dict, else clear it
                last_result = result if isinstance(result, dict) else {}
                logger.debug(f"Updated last_result for step '{step}'.")

        else:
            workflow_func = get_workflow(step)
            if workflow_func:
                workflow_func(current_xyz, config)
            else:
                raise ValueError(f"Step '{step}' is neither a task nor a workflow.")

    logger.info("Workflow complete.")
