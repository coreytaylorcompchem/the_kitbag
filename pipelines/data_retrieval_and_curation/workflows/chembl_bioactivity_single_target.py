from pipeline.workflow_registry import register_workflow
from pipeline.task_registry import get_task

@register_workflow(
    'chembl_bioactivity_single_target',
    description="Retrieve and clean ChEMBL bioactivities and compounds."
)
def run_chembl_bioactivity_workflow(config):
    steps = config.get('workflow', [])
    data = None

    for step in steps:
        print(f"\n>>> Running step: {step}")
        task_func = get_task(step)
        if not task_func:
            raise ValueError(f"Task '{step}' not found in registry.")
        data = task_func(config, data)

    return data
