_task_registry = {}
_task_metadata = {}

def register_task(name: str, category=None, description: str = ""):
    def decorator(func):
        _task_registry[name] = func
        _task_metadata[name] = {
            "description": description,
            "supported_backends": [],
            'category': category or '',
        }
        return func
    return decorator

def get_task(name):
    return _task_registry.get(name)

def list_tasks():
    return list(_task_registry.keys())

def get_task_metadata(name):
    return _task_metadata.get(name, {})

def populate_supported_backends_from_backends(backend_classes: dict):
    """
    Matches each registered task to backends that declare support via 'supported_tasks'.
    This should be called after all tasks and backends have been imported.
    """
    for task_name in _task_registry:
        for backend_name, backend_cls in backend_classes.items():
            supported = getattr(backend_cls, "supported_tasks", [])
            if task_name in supported:
                _task_metadata[task_name]["supported_backends"].append(backend_name)

def finalize_task_registration():
    """
    Should be called after all task modules and backend modules have been loaded.
    """
    from backends import get_backend_classes
    populate_supported_backends_from_backends(get_backend_classes())
