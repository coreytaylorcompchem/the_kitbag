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

def load_all_workflows():
    for _, name, _ in pkgutil.iter_modules(workflows.__path__):
        importlib.import_module(f"{workflows.__name__}.{name}")

# Load all workflows at import time
load_all_workflows()