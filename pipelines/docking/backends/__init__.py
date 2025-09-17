import pkgutil
import importlib
from pathlib import Path

backend_registry = {}

def discover_backends():
    global backend_registry
    backend_registry = {}

    # The directory of this __init__.py file
    backends_path = Path(__file__).parent

    # The package name for importlib
    package_name = __name__  # 'backends'

    # Iterate over all modules in the 'backends' directory
    for finder, name, ispkg in pkgutil.iter_modules([backends_path]):
        # Skip __init__.py or any non-backend files if you want (optional)
        if name.startswith("_"):
            continue
        
        module_name = f"{package_name}.{name}"
        module = importlib.import_module(module_name)

        # Look for attribute 'backend_class' first
        if hasattr(module, "backend_class"):
            cls = getattr(module, "backend_class")
            key = cls.__name__.lower().replace("backend", "")
            backend_registry[key] = cls
            continue

        # Otherwise, find any class ending with 'Backend'
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

def backends_supporting_task(task_name):
    supported = []
    for name, backend_class in backend_registry.items():
        if hasattr(backend_class, task_name):
            supported.append(name)
    return supported

# Discover backends immediately when module is loaded
discover_backends()
