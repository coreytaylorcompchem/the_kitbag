from openmm.app import *
from openmm import *
from openmm.unit import *
import os
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

from backends.utils.add_terminal_caps import CapTermini

class OpenMMBackend:
    def __init__(self, config):
        self.config = config
        self.platform = Platform.getPlatformByName(config["simulation"]["platform"])

    def load_pdb(self, pdb_file):
        return PDBFile(pdb_file)

    def create_system(self, modeller, forcefield_files):
        ff = ForceField(*forcefield_files)
        return ff.createSystem(modeller.topology, nonbondedMethod=PME, constraints=HBonds)

    def add_solvent(self, modeller, ionic_strength, box_padding):
        forcefield = ForceField(*self.config["system"]["forcefield"])
        modeller.addSolvent(forcefield, model='tip3p', ionicStrength=ionic_strength*molar, padding=box_padding*nanometers)
        return modeller
    
    def fix_pdb(self, pdb_file, pH=7.0):
        fixer = PDBFixer(filename=pdb_file)

        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        
        missing_residues = fixer.missingResidues
        missing_atoms = {res: [atom.name for atom in atoms] for res, atoms in fixer.missingAtoms.items()}
        missing_termini = fixer.missingTerminals

        logger.info("Before fixing:")
        logger.info(f"Missing residues: {missing_residues}")
        logger.info(f"Missing atoms: {missing_atoms}")
        logger.info(f"Missing terminals: {missing_termini}")

        logger.info(f"Fixing missing atoms and adding Hydrogens")

        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=pH)

        logger.info(f"Done fixing missing atoms and adding Hydrogens")

        # print("After fixing:")
        # for atom in fixer.topology.atoms():
        #     if atom.name in {"CG", "CD", "OE1", "OE2", "ND2", "OD1"}:
        #         print(f"Found atom {atom.name} in residue {atom.residue.index} ({atom.residue.name})")

        fixed_pdb_path = pdb_file.replace(".pdb", "_fixed.pdb")
        with open(fixed_pdb_path, 'w') as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)

        return fixed_pdb_path

    def mutate_residue_in_pdb(self, input_pdb_path, output_pdb_path, from_resname, to_resname):
        """
        Helper to convert specified residues (e.g., SEP → SER) and removes non-standard atoms (e.g., phosphate group).
        TODO: make this work with other non-standard residues.
        It'll take in other residues than SEP/SER but it won't fix them correctly at present.
        """
        phosphate_atoms = {"P", "OP1", "OP2", "OP3", "O1P", "O2P", "O3P"}

        with open(input_pdb_path, 'r') as f:
            lines = f.readlines()

        new_lines = []
        for line in lines:
            if line.startswith(("ATOM", "HETATM")):
                resname = line[17:20]
                atom_name = line[12:16].strip()

                if resname == from_resname:
                    if atom_name in phosphate_atoms:
                        continue  # Skip phosphate atoms
                    # Simple replacement of residue name (e.g., SEP → SER)
                    line = line[:17] + to_resname.ljust(3) + line[20:]

            new_lines.append(line)

        with open(output_pdb_path, 'w') as f:
            f.writelines(new_lines)

        return output_pdb_path
    
    def prepare_system(self):
        # Mutate SEP -> SER
        mutated_pdb = self.config["system"]["pdb_file"].replace(".pdb", "_mutated.pdb")
        self.mutate_residue_in_pdb(self.config["system"]["pdb_file"], mutated_pdb, 'SEP', 'SER')

        # PDBFixer run
        fixed_pdb_file = self.fix_pdb(mutated_pdb, pH=self.config["simulation"]["pH"])
        pdb = PDBFile(fixed_pdb_file)

        # Delete anything other than protein and HOH.
        atoms_to_keep = [atom for atom in pdb.topology.atoms()
                 if (atom.residue.name in ('HOH', 'WAT'))]

        modeller = Modeller(pdb.topology, pdb.positions)
        modeller.delete([atom for atom in modeller.topology.atoms() if atom not in atoms_to_keep])
        
        # Solvate system
        forcefield_files = self.config["system"]["forcefield"]
        forcefield = ForceField(*forcefield_files)

        ionic_strength = self.config["system"].get("ionic_strength", 0.0)
        box_padding = self.config["system"].get("box_padding", 1.0)

        logger.info(f"Solvating the system for heating/equilibration.")

        modeller = self.add_solvent(modeller, ionic_strength=ionic_strength, box_padding=box_padding)

        logger.info(f"Done solvating the system for heating/equilibration.")

        # Compile the complete system
        self.system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=PME,
            constraints=HBonds
        )

        self.topology = modeller.topology
        self.positions = modeller.positions

