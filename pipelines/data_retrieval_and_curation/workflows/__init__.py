import importlib
import pkgutil

# Registry for workflows
_workflow_registry = {}

def register_workflow(name, description=None):
    def decorator(func):
        _workflow_registry[name.lower()] = {"func": func, "description": description}
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
    import workflows  # make sure this is referencing the module, not just the folder
    for _, name, _ in pkgutil.iter_modules(workflows.__path__):
        try:
            importlib.import_module(f"{workflows.__name__}.{name}")
        except Exception as e:
            print(f"[load_all_workflows] ⚠️ Failed to import workflow module '{name}': {e}")