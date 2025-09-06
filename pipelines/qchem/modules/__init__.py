import importlib
import pkgutil
import modules

# Auto-load all task modules on import
def load_all_tasks():
    for _, name, _ in pkgutil.iter_modules(modules.__path__):
        importlib.import_module(f"{modules.__name__}.{name}")

load_all_tasks()
