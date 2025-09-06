from abc import ABC, abstractmethod

class BaseBackend(ABC):
    def __init__(self, config):
        self.config = config

    @abstractmethod
    def run(self, xyz_file: str, task: str):
        pass