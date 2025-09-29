from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path

def load_peptide_sequences(sequence_file: Path):
    """
    Load peptide sequences from a file. Supports FASTA or plain text formats.
    Returns: list of SeqRecord
    """
    sequences = []

    # Try FASTA parsing first
    try:
        records = list(SeqIO.parse(sequence_file, "fasta"))
        if records:
            return records
    except Exception:
        pass  # fallback to plain text

    # Plain text fallback
    with open(sequence_file, "r") as f:
        lines = [line.strip() for line in f if line.strip()]
        for idx, seq in enumerate(lines):
            record = SeqRecord(Seq(seq), id=f"peptide_{idx+1}")
            sequences.append(record)

    return sequences
