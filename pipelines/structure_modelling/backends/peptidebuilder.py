import os
from pathlib import Path
import logging
import yaml
import subprocess

import PeptideBuilder
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import EditableMol
from Bio.PDB import PDBIO
from Bio import SeqIO

logger = logging.getLogger(__name__)

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdchem

from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit import Chem
from rdkit.Chem import AllChem

from openmm.app import PDBFile, Modeller, ForceField
from openmm.unit import nanometer
from pdbfixer import PDBFixer

def add_caps(peptide_pdb_path, output_pdb_path):
    # Load the peptide using PDBFixer
    fixer = PDBFixer(filename=peptide_pdb_path)

    # Try to detect and fix missing parts
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    # Add hydrogens and use neutral terminal patches
    fixer.addMissingHydrogens(pH=7.0)
    fixer.addMissingTerminals = True  # Enable neutral capping

    # Convert to OpenMM objects
    modeller = Modeller(fixer.topology, fixer.positions)

    # Use a standard force field for compatibility
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')
    modeller.addHydrogens(forcefield)  # Ensure Hs are consistent with FF

    # Save final capped structure
    with open(output_pdb_path, 'w') as out_file:
        PDBFile.writeFile(modeller.topology, modeller.positions, out_file)

SECONDARY_STRUCTURES = {
    'alpha_helix': {'phi': -57.0, 'psi': -47.0},
    'beta_sheet': {'phi': -139.0, 'psi': 135.0},
    'polyproline_II': {'phi': -75.0, 'psi': 145.0},
    'left_alpha_helix': {'phi': 57.0, 'psi': 47.0},
    'coil': {'phi': -90.0, 'psi': 135.0},
}

class PeptideBuilderBackend:
    def __init__(self, sequence_file=None, output_dir=None, cap_termini=True):
        self.sequence_file = Path(sequence_file) if sequence_file else None
        self.output_dir = Path(output_dir or ".")
        self.cap_termini = cap_termini
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def read_fasta_like(self):
        """
        Reads your >pep1 style FASTA and yields (name, sequence)
        """
        records = list(SeqIO.parse(str(self.sequence_file), "fasta"))
        for rec in records:
            yield rec.id, str(rec.seq)

    def build_peptide(self, seq, out_pdb, structure_type=None):
        angles = SECONDARY_STRUCTURES.get(structure_type, SECONDARY_STRUCTURES['coil'])
        phi = [angles['phi']] * len(seq)
        psi = [angles['psi']] * len(seq)
        psi_im1 = [psi[0]] + psi[:-1]

        structure = PeptideBuilder.make_structure(seq, phi, psi_im1)

        io = PDBIO()
        io.set_structure(structure)
        io.save(str(out_pdb))
        logger.info(f"Built peptide ({structure_type}) saved: {out_pdb}")

    def run(self):
        results = {}
        for name, seq in self.read_fasta_like():
            pep_pdb = self.output_dir / f"{name}_built.pdb"
            capped_pdb = self.output_dir / f"{name}_capped.pdb"
            
            self.build_peptide(seq, pep_pdb)

            if self.cap_termini:
                add_caps(str(pep_pdb), str(capped_pdb))
                final_pdb = capped_pdb
            else:
                final_pdb = pep_pdb

            results[name] = {
                "capped_pdb": str(final_pdb),
            }

        return results