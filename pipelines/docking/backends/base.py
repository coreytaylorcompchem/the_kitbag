from abc import ABC, abstractmethod
from pathlib import Path

class BaseBackend(ABC):
    def __init__(self, executable_path: str, use_gpu: bool = False):
        self.executable_path = Path(executable_path)
        self.use_gpu = use_gpu
        self.cache = {}

    @abstractmethod
    def dock(self, ligand: dict, config: dict):
        """Perform docking for a given ligand and config."""
        pass

    def prepare(self):
        """Optional preparation step (e.g., for the executable)."""
        pass

    def supports_gpu(self) -> bool:
        """Return True if this backend supports GPU acceleration."""
        return self.use_gpu
