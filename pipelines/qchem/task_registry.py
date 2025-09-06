# Register all qchem tasks and allows us to store metadata about the tasks

from backends import backends_supporting_task

_task_registry = {}

def register_task(name, description=None, modifies_geometry=False, supported_backends=None):
    # If not given, detect which backends implement it
    if supported_backends is None:
        supported_backends = backends_supporting_task(name)

    def wrapper(func):
        _task_registry[name.lower()] = {
            'func': func,
            'description': description or '',
            'modifies_geometry': modifies_geometry,
            'supported_backends': supported_backends,
        }
        return func
    return wrapper

used_task_names = set()

def get_task(name):
    name = name.lower()
    used_task_names.add(name)
    entry = _task_registry.get(name)
    return entry['func'] if entry else None

def list_tasks():
    return list(_task_registry.keys())

def get_task_metadata(name):
    return _task_registry.get(name.lower(), {})

def get_unregistered_tasks_used():
    return [name for name in used_task_names if name not in _task_registry]
