import sys
from pathlib import Path
import modeller
import modeller.automodel
from Bio import SeqIO, pairwise2
from Bio.PDB import PDBParser, PPBuilder

def get_first_residue_and_chain(pdb_file):
    """Get the first residue number and chain ID from the PDB file."""
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM"):
                chain_id = line[21].strip() or 'A'
                res_num = int(line[22:26].strip())
                return res_num, chain_id
    return 1, 'A'

def extract_pdb_sequence(pdb_file, chain_id='A'):
    """Extract amino acid sequence from specified chain in PDB file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", pdb_file)
    model = structure[0]
    for chain in model:
        if chain.id == chain_id:
            ppb = PPBuilder()
            peptides = ppb.build_peptides(chain)
            sequence = "".join(str(pp.get_sequence()) for pp in peptides)
            return sequence
    raise ValueError(f"Chain {chain_id} not found in PDB")

def load_fasta_sequence(fasta_file):
    """Load sequence from a FASTA file."""
    record = SeqIO.read(fasta_file, "fasta")
    return str(record.seq)

def align_sequences(pdb_seq, fasta_seq):
    """Perform global pairwise alignment between PDB and FASTA sequences."""
    alignments = pairwise2.align.globalxx(pdb_seq, fasta_seq)
    # Take the best alignment
    aln_pdb, aln_fasta, score, begin, end = alignments[0]
    return aln_pdb, aln_fasta

def write_pir_alignment(pir_path, code, pdb_aln, fasta_aln, start_res=1, chain_id='A'):
    """Write PIR alignment file for Modeller using aligned sequences."""
    with open(pir_path, 'w') as f:
        # Structure block
        f.write(f">P1;{code}\n")
        f.write(f"structureX:{code}:{start_res}:{chain_id}:{start_res + len(pdb_aln.replace('-', '')) - 1}:{chain_id}:::0.00: 0.00\n")
        f.write(pdb_aln + "*\n")
        # Sequence block
        f.write(f">P1;{code}_seq\n")
        f.write(f"sequence:{code}_seq:::::::0.00: 0.00\n")
        f.write(fasta_aln + "*\n")

def run_modeller(pdb_file, fasta_file):
    """Run Modeller to rebuild missing atoms using known structure and sequence."""
    env = modeller.Environ()
    env.io.atom_files_directory.append(str(Path(pdb_file).parent.resolve()))

    pdb_path = Path(pdb_file)
    pdb_stem = pdb_path.stem
    ali_path = f"{pdb_stem}.ali"

    start_res, chain_id = get_first_residue_and_chain(pdb_file)
    print(f"Using start residue {start_res} and chain ID '{chain_id}'")

    pdb_seq = extract_pdb_sequence(pdb_file, chain_id)
    fasta_seq = load_fasta_sequence(fasta_file)

    pdb_aln, fasta_aln = align_sequences(pdb_seq, fasta_seq)
    write_pir_alignment(ali_path, pdb_stem, pdb_aln, fasta_aln, start_res, chain_id)

    a = modeller.automodel.AutoModel(env,
                                     alnfile=ali_path,
                                     knowns=pdb_stem,
                                     sequence=f"{pdb_stem}_seq")
    a.starting_model = 1
    a.ending_model = 1

    print("Running Modeller to build missing atoms...")
    a.make()

    print("Modelling complete.")
    output_model = f"{pdb_stem}_seq.B99990001.pdb"
    print(f"Rebuilt structure written to: {output_model}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python modeller_test.py your_pdb_file.pdb your_sequence.fasta")
        sys.exit(1)

    run_modeller(sys.argv[1], sys.argv[2])

