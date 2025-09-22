from backends.psi4 import Psi4Backend
from backends.xtb import XTBBackend
# from backends.orca import OrcaBackend

def get_backend_class(name):
    name = name.lower()
    if name == "psi4":
        return Psi4Backend
    elif name == "xtb":
        return XTBBackend
    # elif name == "orca":
    #     return OrcaBackend
    else:
        raise ValueError(f"Unknown backend: {name}")
