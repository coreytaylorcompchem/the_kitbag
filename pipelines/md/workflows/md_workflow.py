from workflows import register_workflow
from backends.openmm_backend import OpenMMBackend
from modules.md import MDWorkflow

@register_workflow(
    "molecular_dynamics",
    description="Perform MD."
)
def run_basic_md_workflow(config):
    backend = OpenMMBackend(config)
    workflow = MDWorkflow(backend, config)

    workflow.prepare_system()
    workflow.setup_simulation()
    workflow.minimize()
    workflow.heat_and_equilibrate()
    workflow.run_production()
