from backends.utils.add_terminal_caps import CapTermini
from backends.utils.load_ligand import load_ligand_from_sdf

from openmm.app import *
from openmm import *
from openmm.unit import *
import os
from pdbfixer import PDBFixer
from openmm.app import PDBFile, Modeller, ForceField
from pipeline.logger import setup_logger

from simtk.openmm.app import PDBFile

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

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
        config = self.config
        ligand_file = config["system"].get("ligand_file")

        # Mutate SEP -> SER
        mutated_pdb = config["system"]["pdb_file"].replace(".pdb", "_mutated.pdb")
        self.mutate_residue_in_pdb(config["system"]["pdb_file"], mutated_pdb, 'SEP', 'SER')

        # Fix with PDBFixer
        fixed_pdb_file = self.fix_pdb(mutated_pdb, pH=config["simulation"]["pH"])
        pdb = PDBFile(fixed_pdb_file)

        modeller = Modeller(pdb.topology, pdb.positions)

        # Remove non-protein, non-water atoms (e.g., crystallographic ligands, ions)
        atoms_to_keep = [
            atom for atom in modeller.topology.atoms()
            if atom.residue.name in ('HOH', 'WAT') or atom.residue.chain.id in {'A'}
        ]
        modeller.delete([atom for atom in modeller.topology.atoms() if atom not in atoms_to_keep])
        logger.info("Removing non-protein atoms, retaining any crystal waters.")

        # If ligand is provided in the yaml, load and add it to the system
        if ligand_file:
            logger.info(f"Ligand file specified: {ligand_file}")
            ligand_pdb_path = load_ligand_from_sdf(ligand_file)
            ligand = PDBFile(ligand_pdb_path)
            modeller.add(ligand.topology, ligand.positions)
            logger.info("Ligand merged into system.")

        # Solvate
        forcefield_files = config["system"]["forcefield"]
        ionic_strength = config["system"].get("ionic_strength", 0.0)
        box_padding = config["system"].get("box_padding", 1.0)

        modeller = self.add_solvent(modeller, ionic_strength=ionic_strength, box_padding=box_padding)

        # Final system creation
        forcefield = ForceField(*forcefield_files)
        self.system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=PME,
            constraints=HBonds
        )

        self.topology = modeller.topology
        self.positions = modeller.positions

        # Determine output dir to save topology
        output_trajectory = self.config["production"]["output_trajectory"]
        output_dir = os.path.dirname(output_trajectory)
        
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        topology_filename = "topology.pdb"
        topology_path = os.path.join(output_dir, topology_filename)

        # Save topology
        with open(topology_path, "w") as f:
            PDBFile.writeFile(self.topology, self.positions, f)

        logger.info(f"Saved prepared system topology to {topology_path}")


