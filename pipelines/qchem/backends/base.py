from abc import ABC

class BaseBackend(ABC):
    def __init__(self, config):
        self.config = config

    def run(self, xyz_file: str, task: str, task_config: dict = None):
        """
        Optional unified entry point to run a backend task by name.
        Dispatches to a method if one exists (e.g., 'optimise', 'single_point').
        """
        task_method = getattr(self, task, None)
        if task_method is None or not callable(task_method):
            raise NotImplementedError(f"Backend '{type(self).__name__}' does not support task: '{task}'")
        return task_method(xyz_file, task_config or {})
