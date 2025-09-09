_WORKFLOW_REGISTRY = {}

def register_workflow(name, description=""):
    def decorator(func):
        _WORKFLOW_REGISTRY[name] = func
        return func
    return decorator

def get_workflow(name):
    return _WORKFLOW_REGISTRY.get(name)
