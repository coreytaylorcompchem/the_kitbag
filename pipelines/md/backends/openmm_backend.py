from openmm.app import *
from openmm import *
from openmm.unit import *
import os
from pdbfixer import PDBFixer
from openmm.app import PDBFile

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
        return modeller  # Return the updated modeller object
    
    def fix_pdb(self, pdb_file, pH=7.0):
        fixer = PDBFixer(filename=pdb_file)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=pH)
        print("Missing terminals:", fixer.missingTerminals)
        print("Missing residues:", fixer.missingResidues)
        print("Missing atoms:", fixer.missingAtoms)

        # if self.config["system"].get("add_terminal_caps", True):
        #     fixer.addMissingTerminals()
        
        fixed_pdb_path = pdb_file.replace(".pdb", "_fixed.pdb")
        with open(fixed_pdb_path, 'w') as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)
        return fixed_pdb_path

    def mutate_residue_in_pdb(self, input_pdb_path, output_pdb_path, from_resname, to_resname):
        with open(input_pdb_path, 'r') as f:
            lines = f.readlines()

        new_lines = []
        for line in lines:
            if line.startswith(("ATOM", "HETATM")):
                resname = line[17:20]
                if resname == from_resname:
                    line = line[:17] + to_resname + line[20:]
            new_lines.append(line)

        with open(output_pdb_path, 'w') as f:
            f.writelines(new_lines)

        return output_pdb_path
    
    def prepare_system(self):
        # Fix non=standard residues
        mutated_pdb = self.config["system"]["pdb_file"].replace(".pdb", "_mutated.pdb")
        self.mutate_residue_in_pdb(self.config["system"]["pdb_file"], mutated_pdb, 'SEP', 'SER')

        # Step 2: Cap termini
        if self.config["system"].get("add_terminal_caps", True):
            capper = CapTermini(mutated_pdb)
            capped_pdb = capper.cap()
        else:
            capped_pdb = mutated_pdb
        
        fixed_pdb_file = self.fix_pdb(mutated_pdb, pH=7.0)
        # First, fix the input PDB to add missing residues, atoms, and hydrogens, and protonate
        fixed_pdb_file = self.fix_pdb(self.config["system"]["pdb_file"], pH=7.0)

        # Now load the fixed PDB, which should have all hydrogens and correct protonation
        pdb = PDBFile(fixed_pdb_file)
        modeller = Modeller(pdb.topology, pdb.positions)

        forcefield_files = self.config["system"]["forcefield"]
        forcefield = ForceField(*forcefield_files)

        # Add solvent box if specified in config (use defaults if missing)
        ionic_strength = self.config["system"].get("ionic_strength", 0.0)
        box_padding = self.config["system"].get("box_padding", 1.0)
        modeller = self.add_solvent(modeller, ionic_strength=ionic_strength, box_padding=box_padding)

        # Create system with constraints and PME for long-range electrostatics
        self.system = forcefield.createSystem(modeller.topology,
                                            nonbondedMethod=PME,
                                            constraints=HBonds)

        # Store updated topology and positions for later use
        self.topology = modeller.topology
        self.positions = modeller.positions
