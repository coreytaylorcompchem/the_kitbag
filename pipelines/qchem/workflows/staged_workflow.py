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
    last_result = None

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

            # Pass wavefunction object if available and step supports it
            if last_result is not None and step == "mesp":
                input_for_task = last_result.get("wfn", current_xyz)

                # If method/basis not set in mesp config, pull from previous step
                if "method" not in step_config and "method" in last_result:
                    step_config["method"] = last_result["method"]
                if "basis" not in step_config and "basis" in last_result:
                    step_config["basis"] = last_result["basis"]
            else:
                input_for_task = current_xyz

            result = task_func(backend, input_for_task, step_config, config)

            if step in ('single_point', 'optimise'):
                # Expecting dict: {'energy': ..., 'wfn': ..., 'method': ..., 'basis': ...}
                if isinstance(result, dict):
                    last_result = result
                    if "wfn" in result and result["wfn"]:
                        current_xyz = result["wfn"].molecule().save_string_xyz()
                else:
                    last_result = None

            elif getattr(task_func, 'modifies_geometry', False) and isinstance(result, str):
                current_xyz = result
                last_result = {"geometry": result}
            else:
                last_result = {}

        else:
            workflow_func = get_workflow(step)
            if workflow_func:
                workflow_func(current_xyz, config)
            else:
                raise ValueError(f"Step '{step}' is neither a task nor a workflow.")

    logger.info("Workflow complete.")
