from task_registry import register_task

@register_task(
    'torsion_scan',
    description="Do a torsion scan of a bond.",
    modifies_geometry=False,
)
def run(backend, xyz_file, config):
    print("Running torsion scan...")
    # ADD THE REST