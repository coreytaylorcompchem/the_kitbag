import builtins

import os
import sys
import subprocess
import threading
import time
import itertools

import numpy as np

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

from backends.utils.orient_gpcr import align_to_opm_reference, orient_gpcr_with_ligand
from backends.utils.transforms import compute_rigid_transform

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

class Spinner:
    def __init__(self, message="Working"):
        self.message = message
        self._stop = threading.Event()

    def start(self):
        def run():
            for c in itertools.cycle("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"):
                if self._stop.is_set():
                    break
                sys.stdout.write(f"\r{self.message} {c}")
                sys.stdout.flush()
                time.sleep(0.1)
            sys.stdout.write("\r" + " " * (len(self.message) + 2) + "\r")

        self.thread = threading.Thread(target=run)
        self.thread.start()

    def stop(self):
        self._stop.set()
        self.thread.join()

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

    def cap_internal_chain_breaks(self, modeller, dist_threshold=0.16*unit.nanometer):
        """
        Split chains at internal breaks, add OXT to upstream C-terminal,
        and rebuild downstream residues with N-terminal Hs.
        Also ensures true chain termini get proper OXT.
        """
        from openmm.app import Topology, Element, Modeller
        from openmm import Vec3, unit
        import numpy as np
        import string

        logger.info("Checking for internal protein chain breaks and fixing termini...")

        old_top = modeller.topology
        old_pos = modeller.positions
        pos_np = [p.value_in_unit(unit.nanometer) for p in old_pos]

        breakpoints = []

        # Detect internal chain breaks
        for chain in old_top.chains():
            residues = list(chain.residues())
            for i in range(len(residues) - 1):
                res_curr = residues[i]
                res_next = residues[i + 1]
                c_atom = next((a for a in res_curr.atoms() if a.name == "C"), None)
                n_atom = next((a for a in res_next.atoms() if a.name == "N"), None)
                if c_atom is None or n_atom is None:
                    continue
                dist = np.linalg.norm(pos_np[c_atom.index] - pos_np[n_atom.index])
                if dist > dist_threshold.value_in_unit(unit.nanometer):
                    logger.info(f"Internal chain break detected: {res_curr} -> {res_next} (C–N distance = {dist:.3f} nm)")
                    breakpoints.append((res_curr, res_next))

        used_chain_ids = set(c.id for c in old_top.chains())
        available_chain_ids = [ch for ch in string.ascii_uppercase if ch not in used_chain_ids]

        new_top = Topology()
        new_pos = []
        atom_map = {}

        # Keep track of upstream residues
        upstream_residues = {up for up, _ in breakpoints}

        # Rebuild topology
        for chain in old_top.chains():
            residues = list(chain.residues())
            current_chain = new_top.addChain(chain.id)
            
            for i, res in enumerate(residues):
                # Start new chain if downstream of an internal break
                for upstream, downstream in breakpoints:
                    if res is downstream:
                        if available_chain_ids:
                            new_chain_id = available_chain_ids.pop(0)
                        else:
                            new_chain_id = str(len(new_top.chains()))
                        current_chain = new_top.addChain(new_chain_id)

                new_res = new_top.addResidue(res.name, current_chain, res.id)

                # Determine if this is the true last residue in the chain
                is_chain_terminal = (i == len(residues) - 1)

                # Copy atoms
                for atom in res.atoms():
                    # Skip all HXT atoms entirely
                    if atom.name == "HXT":
                        continue
                    # Skip OXT for internal upstream residues only
                    if res in upstream_residues and atom.name == "OXT":
                        continue
                    # Skip OXT if already handled for terminal residue
                    if is_chain_terminal and atom.name == "OXT":
                        continue

                    new_atom = new_top.addAtom(atom.name, atom.element, new_res)
                    atom_map[atom] = new_atom
                    pos = old_pos[atom.index].value_in_unit(unit.nanometer)
                    new_pos.append(Vec3(*pos))

                # Add OXT for upstream residues (internal breaks)
                if res in upstream_residues:
                    c_atom = next((a for a in res.atoms() if a.name == "C"), None)
                    o_atom = next((a for a in res.atoms() if a.name == "O"), None)
                    if c_atom and o_atom:
                        c = np.array(old_pos[c_atom.index].value_in_unit(unit.nanometer))
                        o = np.array(old_pos[o_atom.index].value_in_unit(unit.nanometer))
                        v_co = o - c
                        v_co /= np.linalg.norm(v_co)
                        bond_length = 0.125
                        oxt_pos = c - v_co * bond_length
                        new_top.addAtom("OXT", Element.getBySymbol("O"), new_res)
                        new_pos.append(Vec3(*oxt_pos))
                        logger.debug(f"Placed OXT for upstream {res.name} {res.id}")

                # Add OXT for true chain terminal if missing
                if is_chain_terminal and not any(a.name == "OXT" for a in new_res.atoms()):
                    c_atom = next((a for a in res.atoms() if a.name == "C"), None)
                    o_atom = next((a for a in res.atoms() if a.name == "O"), None)
                    if c_atom and o_atom:
                        c = np.array(old_pos[c_atom.index].value_in_unit(unit.nanometer))
                        o = np.array(old_pos[o_atom.index].value_in_unit(unit.nanometer))
                        v_co = o - c
                        v_co /= np.linalg.norm(v_co)
                        bond_length = 0.125
                        oxt_pos = c - v_co * bond_length
                        new_top.addAtom("OXT", Element.getBySymbol("O"), new_res)
                        new_pos.append(Vec3(*oxt_pos))
                        logger.debug(f"Placed OXT for terminal {res.name} {res.id}")

        # Copy bonds
        for bond in old_top.bonds():
            a1, a2 = bond
            if a1 in atom_map and a2 in atom_map:
                new_top.addBond(atom_map[a1], atom_map[a2])

        new_modeller = Modeller(new_top, unit.Quantity(new_pos, unit.nanometer))
        new_modeller.addHydrogens()
        logger.info("Chains fixed: OXT added to upstream residues and terminal residues; hydrogens added downstream.")

        return new_modeller

    def _pre_minimise_termini(self, steps=200):
        """
        Pre-minimise terminal caps:
        - Only restrain canonical protein backbone/heavy atoms.
        - N-terminal hydrogens and OXT atoms are free.
        - Works regardless of membrane, ligands, or waters.
        """
        from openmm import CustomExternalForce, VerletIntegrator, Context, LocalEnergyMinimizer
        from openmm import unit, Vec3

        logger.info("Performing short minimisation of termini atoms")

        system = self.system
        topology = self.topology
        positions = self.positions

        # Identify protein residues (canonical AAs)
        canonical_aas = {
            'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
            'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'
        }

        # Collect atoms to restrain (all protein heavy atoms, excluding termini H/OXT)
        restrain_atoms = []
        for chain in topology.chains():
            residues = [r for r in chain.residues() if r.name.upper() in canonical_aas]
            if not residues:
                continue

            # N-terminal: skip hydrogens
            n_term = residues[0]
            n_term_h = {atom.index for atom in n_term.atoms() if atom.element.symbol == 'H'}

            # C-terminal: skip OXT
            c_term = residues[-1]
            c_term_oxt = {atom.index for atom in c_term.atoms() if atom.name == 'OXT'}

            for res in residues:
                for atom in res.atoms():
                    if atom.index not in n_term_h and atom.index not in c_term_oxt:
                        restrain_atoms.append(atom.index)

        logger.info(f"Restraining {len(restrain_atoms)} protein atoms; termini free")

        # Create per-particle harmonic restraint
        restraint = CustomExternalForce("0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
        restraint.addPerParticleParameter("x0")
        restraint.addPerParticleParameter("y0")
        restraint.addPerParticleParameter("z0")
        restraint.addPerParticleParameter("k")

        for atom in topology.atoms():
            pos = positions[atom.index].value_in_unit(unit.nanometer)
            k_val = 5000.0 if atom.index in restrain_atoms else 0.0
            restraint.addParticle(atom.index, [pos[0], pos[1], pos[2], k_val])

        system.addForce(restraint)
        logger.debug("Applied per-particle restraints")

        # Setup integrator and context
        integrator = VerletIntegrator(0.001)  # ps timestep
        context = Context(system, integrator, self.platform)
        context.setPositions(positions)

        # Minimisation in chunks with logging
        log_interval = 5
        remaining_steps = steps
        iteration = 0
        while remaining_steps > 0:
            this_chunk = min(log_interval, remaining_steps)
            LocalEnergyMinimizer.minimize(context, maxIterations=this_chunk)
            iteration += this_chunk
            state = context.getState(getEnergy=True)
            energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            logger.debug(f"Step {iteration}/{steps}: potential energy = {energy:.1f} kJ/mol")
            remaining_steps -= this_chunk

        # Update positions
        state = context.getState(getPositions=True)
        self.positions = state.getPositions()

        # Cleanup
        del context, integrator
        system.removeForce(system.getNumForces() - 1)

        logger.info("Minimisation of termini complete")

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
                    # Replacement residue name (SEP → SER) TODO: make more robustererer
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
        ligand_cfg = cfg.get("ligand", {})
        ligand_resname = ligand_cfg.get("resname")
        ligand_chain = ligand_cfg.get("chain")
        has_ligand = ligand_resname is not None and ligand_chain is not None # guard for no-ligand case
        ligand_file = cfg.get("ligand_file")
        protein_chains = set(cfg.get("protein_chains", []))
        ionic_strength = cfg.get("ionic_strength", 0.0)
        box_padding = cfg.get("box_padding", 1.0)
        ph = cfg.get("pH", 7.4)

        spinner = Spinner("Preparing system")
        spinner.start()

        full_pdb = PDBFile(pdb_file)
        modeller_subset = Modeller(full_pdb.topology, full_pdb.positions)

        residues_to_keep = [
            r for r in modeller_subset.topology.residues()
            if r.chain.id in protein_chains or (
                has_ligand and r.name.upper() == ligand_resname.upper() and r.chain.id == ligand_chain
            )
        ]
        atoms_to_keep = [a for r in residues_to_keep for a in r.atoms()]
        modeller_subset.delete([a for a in modeller_subset.topology.atoms() if a not in atoms_to_keep])

        subset_pdb_path = pdb_file.replace(".pdb", "_subset.pdb")
        with open(subset_pdb_path, "w") as f:
            PDBFile.writeFile(modeller_subset.topology, modeller_subset.positions, f)
        logger.info(f"Saved subset PDB (protein chain(s) + ligand chain) to {subset_pdb_path}")

        # Step 1: Mutate unnatural residues
        mutated_pdb = subset_pdb_path.replace(".pdb", "_mutated.pdb")
        self.mutate_residue_in_pdb(subset_pdb_path, mutated_pdb, 'SEP', 'SER')

        # Step 2: PDBFixer
        fixed_pdb_file = self.fix_pdb(mutated_pdb, pH=ph)
        pdb = PDBFile(fixed_pdb_file)
        modeller = Modeller(pdb.topology, pdb.positions)

        # Cap all internal chain breaks
        modeller = self.cap_internal_chain_breaks(modeller)
        modeller.addHydrogens(pH=ph)

        # Check termini charges
        def check_chain_termini_charges(modeller):
            canonical_aas = {
                'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
                'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'
            }

            chain_charges = {}

            for chain in modeller.topology.chains():
                residues = list(chain.residues())
                if not residues:
                    continue

                # Skip non-protein chains
                if not all(res.name.upper() in canonical_aas for res in residues):
                    continue

                # N-terminal
                n_term = residues[0]
                h_atoms = [atom for atom in n_term.atoms() if atom.element.symbol == 'H']
                n_charge = +1 if len(h_atoms) >= 2 else 0

                # C-terminal
                c_term = residues[-1]
                oxt_atoms = [atom for atom in c_term.atoms() if atom.name == 'OXT']
                c_charge = -1 if oxt_atoms else 0

                chain_charges[chain.id] = n_charge + c_charge

            total_charge = builtins.sum(chain_charges.values())

            logger.debug(f"Chain termini charges (protein only): {chain_charges}")
            logger.debug(f"Total termini contribution to system charge: {total_charge}")

            return total_charge

        net_termini_charge = check_chain_termini_charges(modeller)
        if net_termini_charge != 0:
            logger.warning("System termini are not neutral, add counterions.")
        else:
            logger.info("System termini charges balanced.")

        logger.debug("=== After cap_internal_chain_breaks ===")
        for chain in modeller.topology.chains():
            residues = [res.name for res in chain.residues()]
            logger.debug(f"Chain {chain.id}: {len(residues)} residues -> {residues}")

        # Map YAML protein chains to all fragments created after capping
        protein_chain_fragments = []

        for chain in modeller.topology.chains():
            # If any residue in this chain came from a YAML-specified protein chain...
            if any(res.chain.id in protein_chains for res in chain.residues()):
                protein_chain_fragments.append(chain.id)

        # logger.info(f"Protein chains to keep after capping: {protein_chain_fragments}")

        logger.debug("Residues after PDBFixer:")
        for residue in modeller.topology.residues():
            logger.debug(f"{residue.name} chain {residue.chain.id}")

        # Step 3–4: Ligand detection
        ligand_atoms = [atom for atom in modeller.topology.atoms()
                        if has_ligand and atom.residue.name.upper() == ligand_resname.upper()
                        and atom.residue.chain.id == ligand_chain]

        ligand_atoms = []
        if has_ligand:
            ligand_atoms = [
                atom for atom in modeller.topology.atoms()
                if has_ligand and atom.residue.name.upper() == ligand_resname.upper()
                and atom.residue.chain.id == ligand_chain
            ]

            if not ligand_atoms:
                raise RuntimeError(f"Ligand {ligand_resname} not found in PDB after fixing!")

            logger.info(f"Ligand {ligand_resname} detected in chain {ligand_chain}")
        else:
            logger.info("No ligand specified → running apo workflow")

        logger.info(f"Ligand {ligand_resname} detected in chain {ligand_chain}")

        # Prepare protein+ligand PDB for orientation (fixed)
        # Keep all protein residues, even if chain IDs changed by chain capping
        residues_to_keep = [
            r for r in modeller.topology.residues()
            if not has_ligand or r.name.upper() != ligand_resname.upper()
            or (has_ligand and r.name.upper() == ligand_resname.upper() and r.chain.id == ligand_chain)
            or (has_ligand and r.name.upper() == ligand_resname.upper() and r.chain.id == ligand_chain) 
        ]
        atoms_to_keep = [atom for r in residues_to_keep for atom in r.atoms()]

        stripped_modeller = Modeller(modeller.topology, modeller.positions)
        stripped_modeller.delete([atom for atom in stripped_modeller.topology.atoms() if atom not in atoms_to_keep])

        input_pdb_dir = os.path.dirname(fixed_pdb_file)
        stripped_pdb_path = os.path.join(input_pdb_dir, "protein_plus_ligand.pdb")
        with open(stripped_pdb_path, "w") as f:
            PDBFile.writeFile(stripped_modeller.topology, stripped_modeller.positions, f)

        logger.debug("After creating protein_plus_ligand.pdb")
        for chain in stripped_modeller.topology.chains():
            residues = [res.name for res in chain.residues()]
            logger.debug(f"Chain {chain.id}: {len(residues)} residues -> {residues}")

        # Capture all ligand atoms before orientation

        ligand_coords_pre = None
        if has_ligand:
            ligand_coords_pre = np.array([
                pos.value_in_unit(unit.angstrom)
                for atom, pos in zip(stripped_modeller.topology.atoms(), stripped_modeller.positions)
                if has_ligand and atom.residue.name.upper() == ligand_resname.upper()
            ])

        if has_ligand:
            logger.debug(f"Ligand atoms pre-orientation: {ligand_coords_pre.shape}")
            logger.info(f"Wrote stripped protein+ligand PDB to {stripped_pdb_path}")

        # Step 4: Determine ligand SDF source
        ligand_sdf_path = None

        if has_ligand:
            if ligand_file and os.path.exists(ligand_file):
                ligand_sdf_path = ligand_file
                logger.info(f"Using user-provided ligand SDF: {ligand_sdf_path}")
            else:
                logger.info("No ligand_file provided - extracting ligand from PDB")
                ligand_atom_indices = [atom.index for atom in ligand_atoms]
                ligand_topology = self.extract_sub_topology(modeller.topology, ligand_atom_indices)
                ligand_positions = unit.Quantity(
                    [modeller.positions[i].value_in_unit(unit.nanometer) for i in ligand_atom_indices],
                    unit.nanometer
                )
                ligand_pdb_path = os.path.join(input_pdb_dir, f"{ligand_resname}_extracted.pdb")
                with open(ligand_pdb_path, "w") as f:
                    PDBFile.writeFile(ligand_topology, ligand_positions, f)
                logger.debug(f"Extracted ligand PDB saved to {ligand_pdb_path}")

                ligand_sdf_path = os.path.join(input_pdb_dir, f"{ligand_resname}_extracted.sdf")
                mol = Chem.MolFromPDBFile(ligand_pdb_path, removeHs=False)
                if mol is None:
                    raise RuntimeError("RDKit failed to read ligand from extracted PDB")
                Chem.MolToMolFile(mol, ligand_sdf_path)
                logger.debug(f"Generated ligand SDF: {ligand_sdf_path}")

        # Step 5: Keep only protein + waters
        residues_to_keep = [
            r for r in modeller.topology.residues()
            if r.chain.id in protein_chain_fragments or r.name in ('HOH', 'WAT', 'SOL')
        ]
        atoms_to_keep = [atom for r in residues_to_keep for atom in r.atoms()]
        modeller.delete([atom for atom in modeller.topology.atoms() if atom not in atoms_to_keep])
        logger.debug("Kept protein and waters only; ligand removed for parameterisation")

        logger.debug("=== After removing ligand, keeping protein + waters ===")
        for chain in modeller.topology.chains():
            residues = [res.name for res in chain.residues()]
            logger.debug(f"Chain {chain.id}: {len(residues)} residues -> {residues}")

        ligand_offmol = None
        ligand_positions_oriented = None

        # Step 6: Forcefield
        protein_ff_files = cfg.get("forcefield", ["amber14-all.xml", "amber14/tip3p.xml"])
        ligand_ff_name = cfg.get("ligand_parameters", "openff-2.0.0.offxml")
        if not ligand_ff_name.endswith(".offxml"):
            ligand_ff_name = f"{ligand_ff_name}.offxml"
        forcefield = ForceField(*protein_ff_files)

        # Step 7a: Orient GPCR for membrane

        opm_ref = cfg.get("opm_reference")
        opm_chain = cfg.get("opm_align_chain", list(protein_chains)[0])

        if cfg.get("membrane", False):
            oriented_pdb_path = os.path.join(input_pdb_dir, "protein_plus_ligand_oriented.pdb")
            logger.info("Orienting GPCR for membrane embedding")

            if opm_ref:
                logger.info(f"Using OPM reference for alignment: {opm_ref}")

                align_to_opm_reference(
                    pdb_path=stripped_pdb_path,
                    output_path=oriented_pdb_path,
                    opm_pdb=opm_ref,
                    target_chain=opm_chain,
                    ligand_resnames=[ligand_resname] if has_ligand else None
                )

            else:
                logger.info("No OPM reference provided → falling back to PCA orientation")

                orient_gpcr_with_ligand(
                    pdb_path=stripped_pdb_path,
                    output_path=oriented_pdb_path,
                    ligand_resnames=[ligand_resname] if has_ligand else None
                )

            pdb_oriented = PDBFile(oriented_pdb_path)
            modeller = Modeller(pdb_oriented.topology, pdb_oriented.positions)

            # Detect protein chains post-orient
            protein_chains_after_orient = [
                chain.id for chain in modeller.topology.chains()
                if any(
                        res.name not in ('HOH', 'WAT', 'SOL') and
                        (not has_ligand or res.name.upper() != ligand_resname.upper())
                        for res in chain.residues()
                    )
            ]
            logger.debug(f"Protein chains detected post-orientation: {protein_chains_after_orient}")

            # Keep protein + water + ligand
            residues_to_keep = [
                r for r in modeller.topology.residues()
                if r.chain.id in protein_chains_after_orient
                or r.name in ('HOH', 'WAT', 'SOL')
                or (has_ligand and r.name.upper() == ligand_resname.upper())
            ]
            atoms_to_keep = [a for r in residues_to_keep for a in r.atoms()]
            modeller.delete([a for a in modeller.topology.atoms() if a not in atoms_to_keep])

            # Ligand coordinates post-orientation
            ligand_coords_post = np.array([
                pos.value_in_unit(unit.angstrom)
                for atom, pos in zip(modeller.topology.atoms(), modeller.positions)
                if has_ligand and atom.residue.name.upper() == ligand_resname.upper()
            ])

            R, t = None, None
            if has_ligand:
                R, t = compute_rigid_transform(ligand_coords_pre, ligand_coords_post)
                logger.debug(f"Ligand atoms post-orientation: {ligand_coords_post.shape}")

        # Step 7b: Parameterise ligand
        ligand_offmol = None

        if has_ligand and ligand_sdf_path and os.path.exists(ligand_sdf_path):
            rdkit_supplier = Chem.SDMolSupplier(ligand_sdf_path, removeHs=False)
            rdkit_mol = next((m for m in rdkit_supplier if m is not None), None)
            if rdkit_mol is None:
                raise RuntimeError("RDKit failed to load ligand from SDF")
            ligand_offmol = Molecule.from_rdkit(rdkit_mol, allow_undefined_stereo=True)
            if not ligand_offmol.conformers:
                raise RuntimeError("Ligand SDF has no 3D coordinates")

        # Step 7c: Remove ligand before membrane
        # Use dynamically detected protein chains instead of old IDs
        if has_ligand:
            residues_to_keep = [
                r for r in modeller.topology.residues()
                if r.chain.id in protein_chains_after_orient or r.name in ('HOH', 'WAT', 'SOL')
            ]
            atoms_to_keep = [a for r in residues_to_keep for a in r.atoms()]
            modeller.delete([a for a in modeller.topology.atoms() if a not in atoms_to_keep])
            logger.debug("Removed ligand and other molecules prior to membrane insertion")

        logger.debug(f"Chains retained before membrane insertion: {[c.id for c in modeller.topology.chains()]}")

        logger.debug("Before membrane insertion")
        for chain in modeller.topology.chains():
            residues = [res.name for res in chain.residues()]
            logger.debug(f"Chain {chain.id}: {len(residues)} residues -> {residues}")

        # --- Step 7d: Add membrane (protein-only) ---
        if cfg.get("membrane", False):
            lipid_type = cfg.get("lipid_type", "POPC")
            padding = cfg.get("membrane_padding", 1.0) * unit.nanometer
            logger.info(f"Adding membrane: {lipid_type}")

            # Make a copy for membrane insertion
            temp_modeller = Modeller(modeller.topology, modeller.positions)

            # Keep only canonical protein residues + water ---
            canonical_aas = {
                'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
                'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'
            }
            residues_to_keep = [
                res for res in temp_modeller.topology.residues()
                if res.name.upper() in canonical_aas or res.name in ('HOH', 'WAT', 'SOL')
            ]
            atoms_to_keep = [atom for res in residues_to_keep for atom in res.atoms()]

            # Delete everything else
            temp_modeller.delete([atom for atom in temp_modeller.topology.atoms() if atom not in atoms_to_keep])

            logger.debug(f"Chains retained before membrane insertion: {[c.id for c in temp_modeller.topology.chains()]}")
            logger.debug("=== Before membrane insertion ===")
            for chain in temp_modeller.topology.chains():
                residues = [res.name for res in chain.residues()]
                logger.debug(f"Chain {chain.id}: {len(residues)} residues -> {residues}")

            # Add membrane around protein
            temp_modeller.addMembrane(forcefield, lipidType=lipid_type, minimumPadding=padding)

            # Replace modeller with membrane-inserted protein; ligand will be re-added later
            modeller = temp_modeller
            logger.info("Membrane successfully added around protein")

        logger.debug("=== After membrane insertion and ligand re-addition ===")
        for chain in modeller.topology.chains():
            residues = [res.name for res in chain.residues()]
            logger.debug(f"Chain {chain.id}: {len(residues)} residues -> {residues}")

        # Step 7e: Re-add ligand with full coordinates
        if cfg.get("membrane", False) and has_ligand and ligand_offmol is not None:
            coords_array = ligand_offmol.conformers[0].to('angstrom').magnitude
            coords_array = (R @ coords_array.T).T + t
            ligand_positions_oriented = unit.Quantity([Vec3(*c) for c in coords_array], unit.angstrom).in_units_of(unit.nanometer)
            ligand_topology = ligand_offmol.to_topology().to_openmm()
            modeller.add(ligand_topology, ligand_positions_oriented)
            logger.debug(f"Re-added ligand {ligand_resname} with {len(ligand_positions_oriented)} atoms after membrane insertion")

        # Register ligand template with SMIRNOFF
        if has_ligand and ligand_offmol is not None:
            smirnoff_generator = SMIRNOFFTemplateGenerator(
                molecules=[ligand_offmol],
                forcefield=ligand_ff_name
            )
            forcefield.registerTemplateGenerator(smirnoff_generator.generator)
            logger.debug(f"Registered SMIRNOFF ligand {ligand_resname} with forcefield before solvation")

        if not cfg.get("membrane", False):
            logger.info(f"Adding solvent with padding={box_padding} nm...")
            modeller.addSolvent(
                forcefield,
                model='tip3p',
                ionicStrength=ionic_strength * unit.molar,
                padding=box_padding * unit.nanometer
            )
        else:
            logger.info("Skipping addSolvent because membrane system is already solvated")

        logger.info(f"Solvent addition complete. Total atoms after solvation: {len(list(modeller.topology.atoms()))}")

        logger.info(f"Creating final system ...")

        # Step 9: Create OpenMM system
        self.system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=PME,
            constraints=HBonds
        )
        self.topology = modeller.topology
        self.positions = modeller.positions

        # Short minimisation of terminal caps.
        
        self._pre_minimise_termini()

        logger.info(f"Final system created. Saving topology.")

        # Step 10: Save topology
        output_trajectory = cfg.get("output_trajectory")
        if output_trajectory and not os.path.exists(output_trajectory):
            os.makedirs(output_trajectory)
        topology_path = os.path.join(output_trajectory, "topology.pdb")
        with open(topology_path, "w") as f:
            PDBFile.writeFile(self.topology, self.positions, f, keepIds=True)
        logger.info(f"Saved prepared protein + ligand system topology to {topology_path}")

        from openmm import XmlSerializer

        system_xml = os.path.join(output_trajectory, "system.xml")

        with open(system_xml, "w") as f:
            f.write(XmlSerializer.serialize(self.system))

        logger.info(f"Saved OpenMM system to {system_xml}")

        spinner.stop()

        # Debug; check parameters
        logger.debug("\n=== ForceField contents ===")
        logger.debug("Files loaded: %s", protein_ff_files)
        if ligand_ff_name:
            logger.debug("Ligand FF: %s", ligand_ff_name)

        logger.debug("\nRegistered template generators:")
        for gen in forcefield._templateGenerators:
            logger.debug("  %s", gen)

        logger.debug("\nForces included in ForceField:")
        for force_name in forcefield._forces:
            logger.debug("  %s", force_name)

        logger.debug("\n=== Created System Forces ===")
        for force in self.system.getForces():
            logger.debug("%s", type(force))
            if isinstance(force, openmm.NonbondedForce):
                logger.debug("  Num particles: %d", force.getNumParticles())