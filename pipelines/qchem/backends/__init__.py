from backends.psi4 import Psi4Backend
from backends.xtb import XTBBackend

backend_registry = {
    "psi4": Psi4Backend,
    "xtb": XTBBackend,
    # "orca": ORCABackend
}

def get_backend(name):
    backend_class = backend_registry.get(name.lower())
    if not backend_class:
        raise ValueError(f"Backend '{name}' is not available.")
    return backend_class()

def list_backends():
    return sorted(backend_registry.keys())

def get_backend_class(name):
    return backend_registry.get(name.lower())

def backends_supporting_task(task_name):
    supported = []
    for name, backend_class in backend_registry.items():
        if hasattr(backend_class, task_name):
            supported.append(name)
    return supported