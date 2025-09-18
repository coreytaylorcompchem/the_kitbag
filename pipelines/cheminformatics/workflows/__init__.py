import importlib
import pkgutil
import workflows
from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

# Registry for workflows
_workflow_registry = {}

def register_workflow(name, description=None):
    """
    Decorator to register a workflow function.
    """
    def decorator(func):
        _workflow_registry[name.lower()] = {"func": func, "description": description}
        return func
    return decorator

def get_workflow(name):
    """
    Retrieve a registered workflow function by its name.
    """
    entry = _workflow_registry.get(name.lower())
    return entry["func"] if entry else None

def list_workflows():
    """
    List all registered workflows.
    """
    return list(_workflow_registry.keys())

def get_workflow_metadata(name):
    """
    Retrieve metadata for a registered workflow by name.
    """
    return _workflow_registry.get(name.lower(), {})

def load_all_workflows():
    logger.info("[load_all_workflows] Starting workflow discovery")
    for finder, name, ispkg in pkgutil.walk_packages(workflows.__path__, prefix=workflows.__name__ + "."):
        try:
            importlib.import_module(name)
            logger.info(f"[load_all_workflows] Imported {name}")
        except Exception as e:
            logger.error(f"[load_all_workflows] Failed to import {name}: {e}")
