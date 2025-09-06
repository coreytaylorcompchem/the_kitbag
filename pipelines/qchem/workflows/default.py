from workflows import register_workflow, get_workflow
from task_registry import get_task

@register_workflow(
    'default',
    description="Performs structure optimisation only."
)
def run_default_workflow(backend, xyz_file, config):
    steps = config.get('workflow', [])
    for step in steps:
        print(f"\n>>> Running workflow step: {step}")
        # print("  Checking if task exists...")
        task_func = get_task(step)
        # print(f"  get_task({step!r}) → {task_func}")
        if task_func:
            task_func(backend, xyz_file, config)
        else:
            print("  Task not found. Checking if it's a workflow...")
            workflow_func = get_workflow(step)
            print(f"  get_workflow({step!r}) → {workflow_func}")
            if workflow_func:
                workflow_func(backend, xyz_file, config)
            else:
                raise ValueError(f"Step '{step}' is neither a registered task nor workflow.")



    print("Default workflow complete")