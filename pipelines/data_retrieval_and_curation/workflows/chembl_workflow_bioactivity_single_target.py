from pipeline.workflow_registry import register_workflow
from pipeline.base_workflow import BaseWorkflow

@register_workflow(
    'chembl_bioactivity_single_target',
    description="Retrieve and clean ChEMBL bioactivities and compounds."
)
def run_chembl_bioactivity_workflow(config):
    workflow = BaseWorkflow(config)
    return workflow.run()
