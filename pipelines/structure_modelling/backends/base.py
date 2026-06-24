from pathlib import Path

class BaseBackend:
    def __init__(self, **kwargs):
        self.pose = None
        self.cache = {}

    def load_pdb(self, pdb_path):
        raise NotImplementedError("Subclasses must implement `load_pdb` method.")

    def save_pose(self, path):
        if self.pose:
            self.pose.dump_pdb(path)
        else:
            raise RuntimeError("No pose is loaded in backend.")

class BaseStructureTool:
    """
    Base class for AI-based structure prediction tools.
    Supports multi-chain sequences, multiple seeds, and GPU assignment.
    """

    name: str = "base"

    
    def __init__(self):
        self.sequences = None
        self.cache = {}

        
        """
        Base class for AI-based structure prediction tools.

        Sequences are passed at runtime (run()), not stored at init.
        """

    def run(
        self,
        run_id: int,
        seed: int,
        device: int,
        output_dir: Path,
        refinement: dict | None = None,
    ) -> dict:
        """
        Run the structure prediction.
        Must be implemented by subclasses.

        Parameters
        ----------
        run_id : int
            Identifier for this run (used for caching/output naming)
        seed : int
            Random seed for reproducibility
        device : int
            CUDA device index
        output_dir : Path
            Directory to write outputs
        refinement : dict | None
            Optional refinement parameters (e.g., local structure)
        """
        raise NotImplementedError
