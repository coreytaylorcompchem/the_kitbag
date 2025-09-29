from backends.base import BaseBackend
from pipeline.logger import setup_logger
import pyrosetta

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

class PyRosettaBackend(BaseBackend):
    supported_tasks = [
        "fix_residues",
        "fix_loops",
        "cap_terminals",
        "refine_loops"
    ]

    def __init__(self, rosetta_args: str = "-mute all"):
        super().__init__(executable_path=None)
        pyrosetta.init(extra_options=rosetta_args)
        self.pose = None
        self.sequence = None

    def load_pdb(self, pdb_path):
        from pyrosetta import pose_from_pdb
        logger.info(f"Loading PDB from {pdb_path}")
        self.pose = pose_from_pdb(str(pdb_path))
        return self.pose
