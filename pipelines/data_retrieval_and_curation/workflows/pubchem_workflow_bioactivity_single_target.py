from workflows import register_workflow
from pipeline.base_workflow import BaseWorkflow

@register_workflow(
    'pubchem_bioactivity_single_target',
    description="Retrieve, standardise and collate PubChem bioactivities for a target."
)
def run_pubchem_bioactivity_workflow(config):
    workflow = BaseWorkflow(config)
    return workflow.run()
