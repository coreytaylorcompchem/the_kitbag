import PeptideBuilder
from Bio.PDB import PDBIO

seq = "ACDEFGHIK"

# Default backbone angles (alpha helix)
phi = [-57.0] * len(seq)
psi = [-47.0] * len(seq)

# psi_im1 is psi of previous residue; for first residue, use something like -47 as well
psi_im1 = [psi[0]] + psi[:-1]

structure = PeptideBuilder.make_structure(seq, phi, psi_im1)

io = PDBIO()
io.set_structure(structure)
io.save("peptide.pdb")

