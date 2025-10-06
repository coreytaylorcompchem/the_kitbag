from backends.xtb import XTBBackend

def get_backend_from_config(task_name, config):
    """
    Given a task name and the global config, return the appropriate backend instance.
    """
    task_cfg = config.get(task_name, {})
    backend_name = task_cfg.get("backend")

    if backend_name == "xtb":
        return XTBBackend()
    else:
        raise ValueError(f"Unknown backend '{backend_name}' for task '{task_name}'")