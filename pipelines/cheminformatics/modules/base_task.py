from abc import ABC, abstractmethod

class BaseTask(ABC):
    def __init__(self, config: dict, backend):
        self.config = config
        self.backend = backend

    @abstractmethod
    def run(self, molecules):
        pass
