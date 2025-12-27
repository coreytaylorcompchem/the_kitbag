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

    def __init__(self, sequence: str | dict):
        """
        Parameters
        ----------
        sequence : str or dict
            Either a single-sequence string (legacy) or a dict of chain_id -> sequence.
        """
        if isinstance(sequence, str):
            # Legacy single-sequence behavior
            self.sequences = {"A": sequence}
        elif isinstance(sequence, dict):
            if not all(isinstance(k, str) and isinstance(v, str) for k, v in sequence.items()):
                raise ValueError("Sequences dict must be chain_id -> sequence string")
            self.sequences = sequence
        else:
            raise ValueError("sequence must be a string or a dict of chain_id -> sequence")

    def prepare_input(self, output_dir: Path) -> Path:
        """
        Write the sequences to a multi-chain FASTA file.
        Returns the Path to the FASTA.
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        fasta_file = output_dir / "input.fasta"
        lines = []
        for chain_id, seq in self.sequences.items():
            lines.append(f">{chain_id}")
            lines.append(seq)
        fasta_file.write_text("\n".join(lines))
        return fasta_file

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
