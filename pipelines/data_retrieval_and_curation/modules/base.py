from abc import ABC, abstractmethod

class BaseTask(ABC):
    def __init__(self, config: dict):
        self.config = config

    @abstractmethod
    def run(self, data):
        pass
