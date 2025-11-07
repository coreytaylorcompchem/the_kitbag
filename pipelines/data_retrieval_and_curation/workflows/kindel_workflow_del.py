from workflows import register_workflow
from pipeline.base_workflow import BaseWorkflow

@register_workflow(
    'kindel_del',
    description="Retrieve DNA-Encoded_Library data from Kindel (DDR and MAPK14)."
)
def run_kindel_del_workflow(config):
    workflow = BaseWorkflow(config)
    return workflow.run()
