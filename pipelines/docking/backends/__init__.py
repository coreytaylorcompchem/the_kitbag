from backends.gnina import GninaBackend
# from backends.unidock_backend import UniDockBackend

backend_registry = {
    "gnina": GninaBackend,
    # "unidock": UniDockBackend,
}

def get_backend(name):
    backend_class = backend_registry.get(name.lower())
    if not backend_class:
        raise ValueError(f"Backend '{name}' is not available.")
    return backend_class()

def list_backends():
    return sorted(backend_registry.keys())

def get_backend_class(name):
    return backend_registry.get(name.lower())

def backends_supporting_task(task_name):
    supported = []
    for name, backend_class in backend_registry.items():
        if hasattr(backend_class, task_name):
            supported.append(name)
    return supported

def load_all_backends():
    # Import your backends here to populate the registry
    import backends.gnina  # example
