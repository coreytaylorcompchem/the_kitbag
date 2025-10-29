import pkgutil
import importlib
from pathlib import Path

backend_registry = {}

def discover_backends():
    global backend_registry
    backend_registry = {}

    backends_path = Path(__file__).parent
    package_name = __name__  

    for finder, name, ispkg in pkgutil.iter_modules([backends_path]):
        if name.startswith("_"):
            continue

        module_name = f"{package_name}.{name}"
        module = importlib.import_module(module_name)

        if hasattr(module, "backend_class"):
            cls = getattr(module, "backend_class")
            key = cls.__name__.lower().replace("backend", "")
            backend_registry[key] = cls
            continue

        # Fallback: find any class ending in 'Backend'
        for attr_name in dir(module):
            attr = getattr(module, attr_name)
            if isinstance(attr, type) and attr_name.lower().endswith("backend"):
                key = attr_name.lower().replace("backend", "")
                backend_registry[key] = attr

def get_backend(name, **kwargs):
    backend_class = backend_registry.get(name.lower())
    if not backend_class:
        raise ValueError(f"Backend '{name}' is not available.")
    return backend_class(**kwargs)

def list_backends():
    return sorted(backend_registry.keys())

def get_backend_class(name):
    return backend_registry.get(name.lower())

def get_backend_classes():
    return backend_registry.copy()

def backends_supporting_task(task_name):
    supported = []
    for name, backend_class in backend_registry.items():
        tasks = getattr(backend_class, "supported_tasks", [])
        if task_name in tasks:
            supported.append(name)
    return supported

discover_backends()
