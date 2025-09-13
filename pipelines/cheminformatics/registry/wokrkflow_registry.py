WORKFLOW_REGISTRY = {}

def register_workflow(name):
    def decorator(cls):
        WORKFLOW_REGISTRY[name] = cls
        return cls
    return decorator

def get_workflow(name):
    return WORKFLOW_REGISTRY.get(name)
