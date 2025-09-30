import os
from pathlib import Path
import logging
import yaml

import PeptideBuilder
from Bio.PDB import PDBIO
from Bio import SeqIO

import openmm
from openmm import app, unit, LangevinIntegrator, Platform
from openmm.app import PDBFile, Simulation, Modeller, ForceField
from openmm.app.element import hydrogen, oxygen, nitrogen, carbon
from openmm.app.topology import Topology

logger = logging.getLogger(__name__)

SECONDARY_STRUCTURES = {
    'alpha_helix': {'phi': -57.0, 'psi': -47.0},
    'beta_sheet': {'phi': -139.0, 'psi': 135.0},
    'polyproline_II': {'phi': -75.0, 'psi': 145.0},
    'left_alpha_helix': {'phi': 57.0, 'psi': 47.0},
    'coil': {'phi': -90.0, 'psi': 135.0},
}

class PeptideBuilderBackend:
    def __init__(self, sequence_file=None, output_dir=None, minimization_cfg=None):
        self.sequence_file = Path(sequence_file) if sequence_file else None
        self.output_dir = Path(output_dir or ".")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.minimization_cfg = minimization_cfg or {}

    def read_fasta_like(self):
        """
        Reads your >pep1 style FASTA and yields (name, sequence)
        """
        records = list(SeqIO.parse(str(self.sequence_file), "fasta"))
        for rec in records:
            yield rec.id, str(rec.seq)

    def build_peptide(self, seq, out_pdb, structure_type='alpha_helix'):
        angles = SECONDARY_STRUCTURES.get(structure_type, SECONDARY_STRUCTURES['coil'])
        phi = [angles['phi']] * len(seq)
        psi = [angles['psi']] * len(seq)
        psi_im1 = [psi[0]] + psi[:-1]

        structure = PeptideBuilder.make_structure(seq, phi, psi_im1)

        io = PDBIO()
        io.set_structure(structure)
        io.save(str(out_pdb))
        logger.info(f"Built peptide ({structure_type}) saved: {out_pdb}")

    def minimise_peptide(self, pdb_path, out_pdb):
        """
        Minimizes peptide with OpenMM using standard terminal groups (NH3+, COO-),
        avoiding PDBFixer. This fixes the residue template matching error.
        """
        steps = self.minimization_cfg.get("steps", 5000)
        platform_name = self.minimization_cfg.get("platform", "CPU")
        ff_xml = self.minimization_cfg.get("forcefield", "amber14-all.xml")
        implicit_solvent = self.minimization_cfg.get("implicit_solvent", True)

        # Load heavy-atom-only structure
        pdb = PDBFile(str(pdb_path))

        # --- Rebuild topology explicitly to define termini ---
        new_top = Topology()
        chain = new_top.addChain()
        residue_map = {}
        
        for res in pdb.topology.residues():
            new_res = new_top.addResidue(res.name, chain)
            for atom in res.atoms():
                new_top.addAtom(atom.name, atom.element, new_res)
            residue_map[res] = new_res

        # Reconnect bonds manually
        for bond in pdb.topology.bonds():
            a1, a2 = bond
            new_top.addBond(a1, a2)

        modeller = Modeller(new_top, pdb.positions)

        forcefield = ForceField(ff_xml)

        modeller.addHydrogens(forcefield, pH=7.0)

        if implicit_solvent:
            system = forcefield.createSystem(modeller.topology,
                                            nonbondedMethod=app.NoCutoff,
                                            constraints=app.HBonds,
                                            implicitSolvent=app.GBSAOBCForce)
        else:
            system = forcefield.createSystem(modeller.topology,
                                            nonbondedMethod=app.NoCutoff,
                                            constraints=app.HBonds)

        integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
        platform = Platform.getPlatformByName(platform_name)
        simulation = Simulation(modeller.topology, system, integrator, platform)
        simulation.context.setPositions(modeller.positions)

        logger.info(f"Minimizing {pdb_path.name} with OpenMM (steps={steps})...")
        simulation.minimizeEnergy(tolerance=10, maxIterations=steps)

        state = simulation.context.getState(getPositions=True)
        positions = state.getPositions()

        with open(out_pdb, 'w') as f:
            PDBFile.writeFile(simulation.topology, positions, f)
        logger.info(f"Minimized structure saved: {out_pdb}")


    def run(self):
        results = {}
        for name, seq in self.read_fasta_like():
            pep_pdb = self.output_dir / f"{name}_built.pdb"
            self.build_peptide(seq, pep_pdb)

            min_pdb = self.output_dir / f"{name}_minimized.pdb"
            self.minimise_peptide(pep_pdb, min_pdb)

            results[name] = {
                "built_pdb": str(pep_pdb),
                "minimized_pdb": str(min_pdb),
            }
        return results
