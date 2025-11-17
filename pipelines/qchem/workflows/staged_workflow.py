from workflows import register_workflow, get_workflow
from pipeline.task_registry import get_task
from pipeline.backend_registry import get_backend_class

import os
from pipeline.logger import setup_logger
logger = setup_logger(__name__, debug_mode=False, simple_format=True)


def _write_temp_xyz(xyz_string, step_name, output_dir):
    """
    Writes an XYZ string to a temporary file and returns the file path.
    """
    os.makedirs(output_dir, exist_ok=True)
    temp_path = os.path.join(output_dir, f"{step_name}_intermediate.xyz")
    with open(temp_path, "w") as f:
        f.write(xyz_string)
    logger.debug(f"Wrote intermediate XYZ file: {temp_path}")
    return temp_path


@register_workflow(
    'staged_workflow',
    description="Performs single or multi-stage calculations."
)
def run_staged_workflow(xyz_file, config):

    steps = config.get('workflow', [])
    global_output_cfg = config.get("output", {})
    output_dir = global_output_cfg.get("directory", "results")

    current_xyz_path = xyz_file   # always a FILE path
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
            backend = BackendClass(config)

            # -------- Task input selection --------
            if last_result and step in ("qtaim", "mesp", "basin", "nci"):
                # Wavefunction tasks
                if last_result.get("wfn_file"):
                    input_for_task = last_result["wfn_file"]
                elif last_result.get("wfn"):
                    input_for_task = last_result["wfn"]
                else:
                    input_for_task = current_xyz_path

                # Fill in method/basis if omitted
                for key in ("method", "basis"):
                    if key not in step_config and key in last_result:
                        step_config[key] = last_result[key]

            else:
                # All other tasks (single_point, optimise)
                input_for_task = current_xyz_path

            logger.debug(f"Passing to task '{step}': {type(input_for_task)} - {input_for_task}")

            # -------- Run task --------
            result = task_func(backend, input_for_task, step_config, config)
            logger.debug(f"Result from step '{step}': {result}")

            # -------- Update state --------
            if step in ('single_point', 'optimise'):

                if isinstance(result, dict):
                    last_result = result

                    # 1 — if geometry string available, convert to temp XYZ file
                    if result.get("geometry"):
                        xyz_string = result["geometry"]
                        current_xyz_path = _write_temp_xyz(
                            xyz_string, step, output_dir
                        )

                    # 2 — if only wfn is returned, derive geometry from wfn
                    elif result.get("wfn"):
                        xyz_string = result["wfn"].molecule().save_string_xyz()
                        current_xyz_path = _write_temp_xyz(
                            xyz_string, step, output_dir
                        )

                    # No geometry, keep old
                    else:
                        logger.debug("No geometry or wfn returned; reusing previous XYZ path.")

            # Generic geometry-modifying tasks (rare)
            elif getattr(task_func, 'modifies_geometry', False) and isinstance(result, str):
                current_xyz_path = _write_temp_xyz(result, step, output_dir)
                last_result = {"geometry": result}

            else:
                last_result = result if isinstance(result, dict) else {}

        else:
            workflow_func = get_workflow(step)
            if workflow_func:
                workflow_func(current_xyz_path, config)
            else:
                raise ValueError(f"Step '{step}' is neither a task nor a workflow.")

    logger.info("Workflow complete.")
