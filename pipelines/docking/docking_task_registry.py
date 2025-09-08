_task_registry = {}

def register_task(name, description="", supported_backends=None):
    if supported_backends is None:
        supported_backends = []
    def decorator(func):
        _task_registry[name] = {
            'func': func,
            'description': description,
            'supported_backends': supported_backends
        }
        return func
    return decorator

def get_task(name):
    entry = _task_registry.get(name)
    return entry['func'] if entry else None

def list_tasks():
    return list(_task_registry.keys())

def get_task_metadata(name):
    return _task_registry.get(name, {})

def load_all_tasks():
    # Import all your task modules here to populate _task_registry
    import modules.docking_tasks