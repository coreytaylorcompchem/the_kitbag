import psi4
from task_registry import register_task

@register_task(
    'sapt0',
    description="Do a sapt0 calculation (psi4 only)",
    modifies_geometry=False,
)
def run(backend, xyz_file, config):
    print("Running SAPT0 calculation...")
    
    if backend.__class__.__name__ != 'Psi4Backend':
        raise RuntimeError("SAPT0 is only available in Psi4.")
    mol = psi4.geometry(open(xyz_file).read())
    psi4.energy('sapt0', molecule=mol)