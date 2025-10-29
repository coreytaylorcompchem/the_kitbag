from backends.utils.add_terminal_caps import CapTermini
from backends.utils.load_ligand import load_ligand_from_sdf

import os
import subprocess

from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField as OFFForceField
from openff.toolkit.topology import Topology as OFFTopology
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

from openmm import *
from openmm.unit import *
from openmm import unit as unit 
from openmm.app import Topology as Topology
from openmm.app import PDBFile, Modeller, ForceField, PME, HBonds
from openmm import unit

from simtk.openmm.app import PDBFile

from rdkit import Chem
from rdkit.Chem import AllChem

from pdbfixer import PDBFixer

from backends.utils.orient_gpcr import orient_gpcr_with_ligand

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

class OpenMMBackend:
    supported_tasks = [ 
            "minimize",
            "prepare_system",
            "setup_simulation",
            "heat_and_equilibrate",
            "production",
            "solvent_hbonds"]
    def __init__(self, config):
        self.config = config
        self.platform = Platform.getPlatformByName(config.get("setup_simulation", {}).get("platform", "CPU"))

    def load_pdb(self, pdb_file):
        return PDBFile(pdb_file)

    def create_system(self, modeller, forcefield_files):
        ff = ForceField(*forcefield_files)
        return ff.createSystem(modeller.topology, nonbondedMethod=PME, constraints=HBonds)

    def add_solvent(self, modeller, forcefield=None, ionic_strength=0.0, box_padding=1.0):
        if forcefield is None:
            forcefield = ForceField(*self.config["prepare_system"]["forcefield"])

        modeller.addSolvent(forcefield, model='tip3p', ionicStrength=ionic_strength*unit.molar, padding=box_padding*unit.nanometer)
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

    def extract_sub_topology(self, topology, atom_indices):
        """
        Create a new OpenMM Topology containing only atoms with indices in atom_indices.
        Preserve chains, residues, bonds for those atoms.
        """
        new_topology = Topology()
        atom_map = {}

        # Create chains
        for chain in topology.chains():
            new_chain = new_topology.addChain(chain.id)
            # Add residues and atoms if they have atoms in atom_indices
            for residue in chain.residues():
                residue_atoms = [atom for atom in residue.atoms() if atom.index in atom_indices]
                if residue_atoms:
                    new_residue = new_topology.addResidue(residue.name, new_chain, residue.id)
                    for atom in residue_atoms:
                        new_atom = new_topology.addAtom(atom.name, atom.element, new_residue)
                        atom_map[atom.index] = new_atom

        # Add bonds only if both atoms are in atom_map
        for bond in topology.bonds():
            if bond[0].index in atom_map and bond[1].index in atom_map:
                new_topology.addBond(atom_map[bond[0].index], atom_map[bond[1].index])

        return new_topology

    def prepare_system(self, cfg):

        pdb_file = cfg["pdb_file"]
        forcefield_files = cfg.get("forcefield", ["amber14-all.xml", "amber14/tip3p.xml"])
        ligand_params = cfg.get("ligand_parameters", None)
        ionic_strength = cfg.get("ionic_strength", 0.0)
        box_padding = cfg.get("box_padding", 1.0)
        ph = cfg.get("pH", 7.4)

        # Step 1: Mutate unnatural residues
        # TODO: make this work with others.
        mutated_pdb = pdb_file.replace(".pdb", "_mutated.pdb")
        self.mutate_residue_in_pdb(pdb_file, mutated_pdb, 'SEP', 'SER')

        # Step 2: Fix with PDBFixer
        fixed_pdb_file = self.fix_pdb(mutated_pdb, pH=ph)
        pdb = PDBFile(fixed_pdb_file)
        modeller = Modeller(pdb.topology, pdb.positions)

        # Step 3: Detect ligand residues automatically
        protein_chains = {'A'}  # Only chain A for now
        ligand_residues = [
            residue.name for residue in modeller.topology.residues()
            if residue.name not in ('HOH', 'WAT', 'SOL') and residue.chain.id not in protein_chains
        ]
        ligand_residues = list(set(ligand_residues))
        ligand_sdf_path = None
        ligand_offmol = None

        if ligand_residues:
            ligand_resname = ligand_residues[0]
            logger.info(f"Detected ligand residue name: {ligand_resname}")

            input_pdb_dir = os.path.dirname(fixed_pdb_file)
            stripped_pdb_path = os.path.join(input_pdb_dir, "protein_plus_ligand.pdb")
            ligand_sdf_path = os.path.join(input_pdb_dir, "ligand.sdf")

            # Step 4a: Create stripped PDB (protein + ligand only)
            stripped_modeller = Modeller(pdb.topology, pdb.positions)

            residues_to_keep = set()
            for residue in stripped_modeller.topology.residues():
                if residue.chain.id in protein_chains or residue.name == ligand_resname:
                    residues_to_keep.add(residue)

            atoms_to_keep = [atom for residue in residues_to_keep for atom in residue.atoms()]
            stripped_modeller.delete([atom for atom in stripped_modeller.topology.atoms() if atom not in atoms_to_keep])

            with open(stripped_pdb_path, "w") as f:
                PDBFile.writeFile(stripped_modeller.topology, stripped_modeller.positions, f)
            logger.info(f"Wrote stripped protein + ligand PDB to {stripped_pdb_path}")

            # Step 4b: Use Open Babel to separate and filter small molecules
            # TODO: this needs to be more robust to different PDB structures
            cmd = [
                "obabel", stripped_pdb_path,
                "-O", os.path.join(input_pdb_dir, "ligand_output.sdf"),
                "--separate", "--filter", "atoms<200" #TODO: only works with small molecules
            ]
            logger.debug(f"Running Open Babel command: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                logger.error(f"Open Babel failed:\n{result.stderr}")
                raise RuntimeError("Open Babel ligand extraction failed.")

            # Step 4c: Identify the most likely ligand file
            sdf_candidates = [
                os.path.join(input_pdb_dir, f)
                for f in os.listdir(input_pdb_dir)
                if f.startswith("ligand_output") and f.endswith(".sdf")
            ]

            if not sdf_candidates:
                raise RuntimeError("No SDF files were generated by Open Babel; ligand extraction failed.")

            chosen_sdf = None
            min_atoms = 99999
            for sdf_file in sdf_candidates:
                try:
                    mol = Chem.SDMolSupplier(sdf_file, removeHs=False)[0]
                    if mol is None:
                        continue
                    num_atoms = mol.GetNumAtoms()
                    if 5 < num_atoms < min_atoms:
                        min_atoms = num_atoms
                        chosen_sdf = sdf_file
                except Exception:
                    continue

            if not chosen_sdf:
                raise RuntimeError("Failed to identify ligand SDF among separated files.")

            os.rename(chosen_sdf, ligand_sdf_path)
            for sdf_file in sdf_candidates:
                if os.path.exists(sdf_file) and sdf_file != ligand_sdf_path:
                    os.remove(sdf_file)

            logger.info(f"Ligand extracted and written to {ligand_sdf_path}")

        else:
            logger.warning("No ligand residues detected; proceeding without ligand extraction.")

        #  Step 5: Keep only protein and waters in modeller object
        residues_to_keep = set()
        for residue in modeller.topology.residues():
            if residue.chain.id in protein_chains or residue.name in ('HOH', 'WAT', 'SOL'):
                residues_to_keep.add(residue)

        atoms_to_keep = [atom for residue in residues_to_keep for atom in residue.atoms()]
        modeller.delete([atom for atom in modeller.topology.atoms() if atom not in atoms_to_keep])

        logger.debug("Kept protein and waters only; ligand and other hetero atoms removed.")

        # Step 6: Parameterize ligand and merge
        if ligand_sdf_path and os.path.exists(ligand_sdf_path):
            try:
                rdkit_supplier = Chem.SDMolSupplier(ligand_sdf_path, removeHs=False)
                rdkit_mol = next((m for m in rdkit_supplier if m is not None), None)
                if rdkit_mol is None:
                    raise RuntimeError("RDKit failed to load ligand from SDF")

                ligand_offmol = Molecule.from_rdkit(rdkit_mol, allow_undefined_stereo=True)

                if not ligand_offmol.conformers:
                    raise RuntimeError("Ligand SDF has no 3D coordinates; cannot add to system.")

                # Convert ligand coordinates to OpenMM units
                conf = ligand_offmol.conformers[0]
                coords_array = conf.to('angstrom').magnitude
                ligand_positions = unit.Quantity(coords_array, unit.angstrom).in_units_of(unit.nanometer)

                # Add ligand topology and positions to modeller
                ligand_topology = ligand_offmol.to_topology().to_openmm()
                modeller.add(ligand_topology, ligand_positions)

                logger.info(f"Ligand {ligand_resname} successfully merged into the modeller.")

            except Exception as e:
                logger.exception(f"Failed to load or parameterise ligand: {e}")
                ligand_offmol = None

        # Step 7: Create hybrid Anber/OFF ForceField 
        # Load Amber for protein/water
        logger.debug(f"Calculating forcefield parameters for ligand {ligand_resname}.")

        protein_ff_files = cfg.get("forcefield", ["amber14-all.xml", "amber14/tip3p.xml"])
        ligand_ff_name = cfg.get("ligand_parameters", "openff-2.0.0.offxml")

        # If user wrote "openff-2.0.0" instead of "openff-2.0.0.offxml", fix it
        if not ligand_ff_name.endswith(".offxml"):
            ligand_ff_name = f"{ligand_ff_name}.offxml"

        logger.debug(f"Using protein forcefields: {protein_ff_files}")
        logger.debug(f"Using ligand parameters: {ligand_ff_name}")

        forcefield = ForceField(*protein_ff_files)

        if ligand_offmol is not None:
            smirnoff_generator = SMIRNOFFTemplateGenerator(
                molecules=[ligand_offmol],
                forcefield=ligand_ff_name
            )
            forcefield.registerTemplateGenerator(smirnoff_generator.generator)
            logger.debug("Registered SMIRNOFFTemplateGenerator for ligand.")
        
        # --- Step 7a: Orient protein for membrane
        if cfg.get("membrane", False):
            oriented_pdb_path = os.path.join(input_pdb_dir, "protein_oriented.pdb")
            logger.info("Orienting GPCR for membrane embedding (may take some time)...")
            orient_gpcr_with_ligand(
                pdb_path=stripped_pdb_path,
                output_path=oriented_pdb_path,
                ligand_resnames=[ligand_resname] if ligand_offmol else None
            )

            # Reload oriented structure into modeller
            pdb_oriented = PDBFile(oriented_pdb_path)
            modeller = Modeller(pdb_oriented.topology, pdb_oriented.positions)

            # Remove pre-existing ligand residue from the oriented PDB
            if ligand_offmol:
                atoms_to_delete = [atom for atom in modeller.topology.atoms() if atom.residue.name == ligand_resname]
                if atoms_to_delete:
                    modeller.delete(atoms_to_delete)
                    logger.debug(f"Removed pre-existing ligand {ligand_resname} from oriented PDB")

                # Hack: re-add ligand using OpenFF topology and coordinates
                ligand_topology = ligand_offmol.to_topology().to_openmm()
                conf = ligand_offmol.conformers[0]
                coords_array = conf.to('angstrom').magnitude
                ligand_positions = unit.Quantity(coords_array, unit.angstrom).in_units_of(unit.nanometer)
                modeller.add(ligand_topology, ligand_positions)
                logger.debug(f"Re-added ligand {ligand_resname} using OpenFF template")

        #  Step 7b: Add membrane
        if cfg.get("membrane", False):
            lipid_type = cfg.get("lipid_type", "POPC")
            padding = cfg.get("membrane_padding", 1.0) * unit.nanometer
            logger.info(f"Adding membrane: {lipid_type}")
            modeller.addMembrane(forcefield, lipidType=lipid_type, minimumPadding=padding)

        # Step 8: Solvate system
        ionic_strength = cfg.get("ionic_strength", 0.0)
        box_padding = cfg.get("box_padding", 1.0)
        modeller.addSolvent(
            forcefield,
            model='tip3p',
            ionicStrength=ionic_strength * unit.molar,
            padding=box_padding * unit.nanometer
        )

        # Step 9: Create final OpenMM System
        self.system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=PME,
            constraints=HBonds
        )

        self.topology = modeller.topology
        self.positions = modeller.positions

        # Step 10: Save prepared topology
        output_trajectory = cfg.get("output_trajectory")
        output_dir = os.path.dirname(output_trajectory)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        topology_path = os.path.join(output_dir, "topology.pdb")
        with open(topology_path, "w") as f:
            PDBFile.writeFile(self.topology, self.positions, f, keepIds=True)

        logger.info(f"Saved prepared protein + ligand system topology to {topology_path}")

