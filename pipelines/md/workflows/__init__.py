import os
import importlib.util

_workflow_registry = {}

def register_workflow(name, description=None):
    def decorator(func):
        _workflow_registry[name.lower()] = {"func": func, "description": description}
        return func
    return decorator

def get_workflow(name):
    entry = _workflow_registry.get(name.lower())
    if entry and "func" in entry:
        return entry["func"]
    return None

def list_workflows():
    return {
        name: meta["description"]
        for name, meta in _workflow_registry.items()
    }

def get_workflow_metadata(name):
    return _workflow_registry.get(name.lower(), {})

def load_all_workflows(workflows_dir="workflows"):
    """Recursively load all workflow modules and register them"""
    for root, dirs, files in os.walk(workflows_dir):
        for file in files:
            if file.endswith(".py") and not file.startswith("__"):
                module_path = os.path.join(root, file)
                module_rel_path = os.path.relpath(module_path, workflows_dir)
                module_name = module_rel_path.replace(os.sep, ".").replace(".py", "")
                full_module_name = f"workflows.{module_name}"

                print(f"Importing workflow module: {full_module_name}")  # 👈 Add this

                if full_module_name not in globals():
                    spec = importlib.util.spec_from_file_location(full_module_name, module_path)
                    module = importlib.util.module_from_spec(spec)
                    spec.loader.exec_module(module)

