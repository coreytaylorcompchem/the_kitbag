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
            backend = BackendClass(config)  # global config passed into backend

            # 💡 PASS BOTH step_config and global config
            result = task_func(backend, current_xyz, step_config, config)

            # Chain geometry if needed
            if getattr(task_func, 'modifies_geometry', False) and isinstance(result, str):
                current_xyz = result

        else:
            workflow_func = get_workflow(step)
            if workflow_func:
                workflow_func(current_xyz, config)
            else:
                raise ValueError(f"Step '{step}' is neither a task nor a workflow.")

    logger.info("Workflow complete.")

