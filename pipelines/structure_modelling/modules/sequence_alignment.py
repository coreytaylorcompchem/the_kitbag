from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
from Bio.Align import substitution_matrices

import logging

logger = logging.getLogger(__name__)

def group_missing_residues(missing_indices):
    """
    Group sequential residue indices into continuous segments.
    Example: [10, 11, 12, 20, 21] → [(10, 12), (20, 21)]
    """
    if not missing_indices:
        return []

    groups = []
    start = prev = missing_indices[0]

    for idx in missing_indices[1:]:
        if idx == prev + 1:
            prev = idx
        else:
            groups.append((start, prev))
            start = prev = idx
    groups.append((start, prev))
    return groups

def align_sequences(seq1: str, seq2: str, gap_open: float = -10, gap_extend: float = -0.5):
    """
    Align two protein sequences using global alignment with BLOSUM62.
    Returns the best alignment (seq1_aln, seq2_aln).
    """
    matrix = substitution_matrices.load("BLOSUM62")

    alignments = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)
    best = alignments[0]
    return best.seqA, best.seqB


def get_missing_indices_from_alignment(aligned_pdb_seq, aligned_uniprot_seq):
    """
    Given aligned sequences, return indices in UniProt where residues are missing in the PDB.
    Indexing is 1-based (like PyRosetta/PDB).
    """
    missing_indices = []
    uni_idx = 0
    for i, (pdb_res, uni_res) in enumerate(zip(aligned_pdb_seq, aligned_uniprot_seq)):
        if uni_res != "-":
            uni_idx += 1
            if pdb_res == "-":
                missing_indices.append(uni_idx)
    return missing_indices


def compare_sequences(pdb_seq: str, uniprot_seq: str):
    """
    Align and compare the sequences, returning missing UniProt residue positions.
    """
    logger.info("Aligning PDB sequence to UniProt sequence...")
    aligned_pdb, aligned_uni = align_sequences(pdb_seq, uniprot_seq)

    logger.debug("\n" + format_alignment(aligned_pdb, aligned_uni, 0, 0, len(aligned_pdb)))
    missing = get_missing_indices_from_alignment(aligned_pdb, aligned_uni)
    
    logger.info(f"Found {len(missing)} residues missing in the PDB structure.")
    return missing