from pipeline.task_registry import register_task

@register_task(
    'mesp_map',
    description="Calculates a Molecular Electrostatic Potential Map",
    category='Molecular properties'
)
def run(backend, xyz_file, config):
    print("Calculating molecular electrostatic potential (MESP) map...")
    # generate cube files, call single point, and visualize using backend