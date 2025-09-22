_TASK_REGISTRY = {}

def register_task(name, category=None, description=None, modifies_geometry=True):
    # # If not given, detect which backends implement it
    # if supported_backends is None:
    #     supported_backends = backends_supporting_task(name)

    def wrapper(func):
        _TASK_REGISTRY[name.lower()] = {
            'func': func,
            'category': category or '',
            'description': description or '',
            'modifies_geometry': modifies_geometry,
            # 'supported_backends': supported_backends,
        }
        return func
    return wrapper

def get_task(name):
    """
    Return the registered task function (callable) for a given name.
    """
    entry = _TASK_REGISTRY.get(name.lower())
    if entry is None:
        return None
    return entry["func"] 

def list_tasks():
    """
    List all registered task names.
    """
    return list(_TASK_REGISTRY.keys())

def get_task_metadata(name):
    """
    Get metadata dictionary for a given task name.
    """
    return _TASK_REGISTRY.get(name.lower(), {})

# def get_unregistered_tasks_used():
#     """
#     Return list of task names used but not registered.
#     """
#     return [name for name in used_task_names if name not in _TASK_REGISTRY]