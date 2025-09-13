from abc import ABC, abstractmethod

class BaseWorkflow(ABC):
    def __init__(self, config, backend):
        self.config = config
        self.backend = backend
        self.tasks = []

    @abstractmethod
    def setup(self):
        pass

    def run(self, input_path):
        molecules = self.backend.load_molecules(input_path)
        for task in self.tasks:
            molecules = task.run(molecules)
        return molecules
