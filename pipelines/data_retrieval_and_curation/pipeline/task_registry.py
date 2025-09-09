_TASK_REGISTRY = {}

def register_task(name):
    def wrapper(func):
        _TASK_REGISTRY[name] = func
        return func
    return wrapper

def get_task(name):
    return _TASK_REGISTRY.get(name)
