from workflows import register_workflow
from pipeline.base_workflow import BaseWorkflow

@register_workflow(
    'chembl_bioactivity_single_target',
    description="Retrieve, standardise and collate bioactivities for a single target - CHEMBL."
)
def run_chembl_bioactivity_workflow(config):
    workflow = BaseWorkflow(config)
    return workflow.run()
