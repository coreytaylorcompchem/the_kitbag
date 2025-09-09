from abc import ABC, abstractmethod

class DataBackend(ABC):
    @abstractmethod
    def fetch_bioactivity(self, compound_id: str):
        pass

    @abstractmethod
    def fetch_structure(self, compound_id: str):
        pass
