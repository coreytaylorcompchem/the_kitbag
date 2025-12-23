from pathlib import Path

import openmm
import openmm.app as app
import openmm.unit as unit
from openmm.app import PDBFile, Modeller, Simulation, Element, Topology

def load_protein(pdb_path: Path) -> PDBFile:
    return PDBFile(str(pdb_path))


def add_positional_restraints(system, modeller, atom_indices, k=1000.0):
    """
    Apply harmonic positional restraints to selected atoms.
    Args:
        system: OpenMM System object
        modeller: OpenMM Modeller containing positions
        atom_indices: list of atom indices to restrain
        k: force constant in kJ/mol/nm^2
    """
    # Make sure the positions are Quantities in nm
    positions = modeller.positions
    if not isinstance(positions[0], openmm.unit.Quantity):
        positions = [p * unit.nanometer for p in positions]

    # Create the harmonic restraint force
    force = openmm.CustomExternalForce("0.5*k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    force.addGlobalParameter("k", k * unit.kilojoule_per_mole / unit.nanometer**2)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")

    # Add particles
    for idx in atom_indices:
        pos = positions[idx].value_in_unit(unit.nanometer)
        force.addParticle(idx, [pos[0], pos[1], pos[2]])

    system.addForce(force)