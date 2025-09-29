from backends.base import BaseBackend
from pipeline.logger import setup_logger
from pathlib import Path
import subprocess
import shutil

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

class ColabFoldBackend(BaseBackend):
    supported_tasks = [
        "build_peptide"
    ]

    def __init__(self, **kwargs):
        super().__init__(executable_path="colabfold_batch")  # assumes colabfold_batch is in PATH
        self.num_models = kwargs.get("num_models", 5)
        self.use_templates = kwargs.get("use_templates", False)
        self.amber_relax = kwargs.get("amber_relax", False)
        self.output_dir = None
        self.model_path = None
    
    def build_peptides_from_file(self, sequence_file: Path, output_dir: Path):
        from utils.sequence_loader import load_peptide_sequences

        sequences = load_peptide_sequences(sequence_file)
        if not sequences:
            raise ValueError(f"No sequences found in file: {sequence_file}")

        output_paths = []

        for seq_record in sequences:
            name = seq_record.id or f"peptide_{len(output_paths)}"
            seq = str(seq_record.seq)

            peptide_dir = output_dir / name
            peptide_dir.mkdir(parents=True, exist_ok=True)

            fasta_file = peptide_dir / f"{name}.fasta"
            fasta_file.write_text(f">{name}\n{seq}")

            args = [
                self.executable_path,
                str(fasta_file),
                str(peptide_dir),
            ]
            if not self.use_templates:
                args.append("--no-template")
            if not self.amber_relax:
                args.append("--no-relax")

            logger.info(f"Running ColabFold on: {name}")
            subprocess.run(args, check=True)

            ranked_pdb = peptide_dir / "ranked_0.pdb"
            if not ranked_pdb.exists():
                raise FileNotFoundError(f"No ranked_0.pdb for {name}")
            
            output_paths.append((name, ranked_pdb))

        return output_paths

    def run_colabfold(self, sequence: str, output_dir: Path):
        logger.info("Running ColabFold on peptide sequence...")
        input_file = output_dir / "peptide.fasta"
        input_file.write_text(f">peptide\n{sequence}")

        args = [
            self.executable_path,
            str(input_file),
            str(output_dir),
        ]
        if not self.use_templates:
            args.append("--no-template")
        if not self.amber_relax:
            args.append("--no-relax")

        subprocess.run(args, check=True)

        # Get best-ranked model (ranked_0.pdb)
        ranked_pdb = output_dir / "peptide" / "ranked_0.pdb"
        if not ranked_pdb.exists():
            raise FileNotFoundError("ColabFold did not produce a ranked_0.pdb file.")

        self.model_path = ranked_pdb
        logger.info(f"ColabFold prediction saved to {ranked_pdb}")
        return ranked_pdb

    def build_peptide(self, sequence: str, output_dir: Path):
        self.output_dir = output_dir
        model_path = self.run_colabfold(sequence, output_dir)
        return model_path
