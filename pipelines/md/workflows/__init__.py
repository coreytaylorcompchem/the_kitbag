import importlib
import pkgutil
import workflows

# Registry for workflows
_workflow_registry = {}

def register_workflow(name, description=None):
    def decorator(func):
        _workflow_registry[name.lower()] = {
            "func": func,
            "description": description
        }
        return func
    return decorator

def get_workflow(name):
    entry = _workflow_registry.get(name.lower())
    return entry["func"] if entry else None

def list_workflows():
    return list(_workflow_registry.keys())

def get_workflow_metadata(name):
    return _workflow_registry.get(name.lower(), {})

def load_all_workflows(package=workflows):
    """Recursively load all workflow modules under the given package"""
    for loader, module_name, is_pkg in pkgutil.iter_modules(package.__path__):
        full_module_name = f"{package.__name__}.{module_name}"
        importlib.import_module(full_module_name)
        # Recursively load sub-packages
        if is_pkg:
            subpackage = importlib.import_module(full_module_name)
            load_all_workflows(subpackage)

# Load all workflows at import time
load_all_workflows()
