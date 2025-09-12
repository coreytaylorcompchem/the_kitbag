import importlib
import pkgutil
import tasks

# Auto-load all task modules on import
def load_all_tasks():
    for _, name, _ in pkgutil.iter_modules(tasks.__path__):
        importlib.import_module(f"{tasks.__name__}.{name}")
