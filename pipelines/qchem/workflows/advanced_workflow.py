from workflows import register_workflow
from task_registry import get_task

@register_workflow('advanced', description="optimize → mesp_map → torsion_scan → sapt0")
def run_advanced_workflow(backend, xyz_file, config):
    print("Starting advanced workflow")

    steps = ['optimize', 'mesp_map', 'torsion_scan', 'sapt0']

    for step in steps:
        print(f"\n>>> Running step: {step}")
        task = get_task(step)  
        if not task:
            raise ValueError(f"Step '{step}' is not a registered task.")
        task(backend, xyz_file, config)

    print("Advanced workflow complete")
