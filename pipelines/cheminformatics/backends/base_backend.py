from abc import ABC, abstractmethod

class BaseBackend(ABC):
    @abstractmethod
    def load_molecules(self, source):
        pass

    @abstractmethod
    def compute_descriptors(self, mol):
        pass

    @abstractmethod
    def filter_molecule(self, mol, **kwargs):
        pass
