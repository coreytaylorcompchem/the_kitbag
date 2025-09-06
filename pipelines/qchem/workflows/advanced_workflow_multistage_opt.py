from workflows import register_workflow
from task_registry import get_task
from backends import get_backend

@register_workflow('advanced_multistage_workflow', description="Runs opt → torsion → opt → MESP with different backends.")
def run_custom_workflow(default_backend, xyz_file, config):
    print("[Workflow] Starting custom multi-step workflow")

    steps = config.get('workflow_steps', [])
    if not steps:
        raise ValueError("No workflow_steps defined in config.yaml")

    for idx, step in enumerate(steps):
        task_name = step.get('task')
        backend_name = step.get('backend', config.get('default_backend', 'psi4'))
        params = step.get('params', {})

        print(f"\n>>> Step {idx+1}: {task_name} using {backend_name}")

        task_func = get_task(task_name)
        if not task_func:
            raise ValueError(f"Task '{task_name}' is not registered.")

        backend = get_backend(backend_name)
        merged_config = {**config, **params}

        task_func(backend, xyz_file, merged_config)

    print("[Workflow] Workflow complete.")