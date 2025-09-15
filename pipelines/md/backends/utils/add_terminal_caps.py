from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
import os

class CapTermini:
    def __init__(self, input_pdb):
        self.input_pdb = input_pdb
        self.output_pdb = input_pdb.replace(".pdb", "_capped.pdb")
        self.structure = PDBParser(QUIET=True).get_structure("protein", self.input_pdb)

    def cap(self):
        model = self.structure[0]
        chain = list(model.get_chains())[0]

        residues = list(chain.get_residues())
        n_term = residues[0]
        c_term = residues[-1]

        # Add ACE cap (N-terminal)
        ace = self._create_residue('ACE', n_term.id[1] - 1)
        chain.insert(0, ace)

        # Add NME cap (C-terminal)
        nme = self._create_residue('NME', c_term.id[1] + 1)
        chain.add(nme)

        io = PDBIO()
        io.set_structure(self.structure)
        io.save(self.output_pdb, select=SelectAll())

        return self.output_pdb

    def _create_residue(self, resname, resseq):
        atoms = {
            'ACE': [
                ("C", (0.0, 0.0, 0.0), "C"),
                ("O", (0.0, 0.0, 0.1), "O"),
                ("H1", (0.1, 0.0, 0.0), "H"),
                ("H2", (-0.1, 0.0, 0.0), "H"),
                ("H3", (0.0, 0.1, 0.0), "H")
            ],
            'NME': [
                ("N", (0.0, 0.0, 0.0), "N"),
                ("H", (0.0, 0.1, 0.0), "H"),
                ("C", (0.1, 0.0, 0.0), "C"),
                ("H1", (0.2, 0.0, 0.0), "H"),
                ("H2", (0.1, 0.1, 0.0), "H"),
                ("H3", (0.1, -0.1, 0.0), "H")
            ]
        }

        res = Residue((' ', resseq, ' '), resname, '')

        for atom_name, coord, element in atoms[resname]:
            # 🛡️ Ensure element is valid (1 or 2-letter symbol, uppercase)
            element = element.capitalize()
            if not element or not isinstance(element, str):
                raise ValueError(f"Invalid element '{element}' for atom '{atom_name}'")

            atom = Atom(atom_name, coord, 1.0, 1.0, '', atom_name, 0, element=element)
            res.add(atom)

        return res

class SelectAll(Select):
    def accept_residue(self, residue):
        return True
