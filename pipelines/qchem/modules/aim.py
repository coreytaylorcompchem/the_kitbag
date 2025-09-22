import psi4
from pipeline.task_registry import register_task

@register_task(
    'qtaim',
    description="Do AIM (Atoms In Molecules) analysis.",
    modifies_geometry=False,
    category='EDA'
)
def run(backend, xyz_file, config):
    print("Running QTAIM calculation...")
    
    if backend.__class__.__name__ != 'Psi4Backend':
        raise RuntimeError("SAPT0 is only available in Psi4.")
    mol = psi4.geometry(open(xyz_file).read())
    psi4.energy('sapt0', molecule=mol)