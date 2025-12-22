import os
from pathlib import Path

import openmm
import openmm.app as app
import openmm.unit as unit
from openmm.app import PDBFile, Modeller, Simulation, Element, Topology
from openmmforcefields.generators import SystemGenerator
from openff.toolkit.topology import Molecule
# from openmm.app.element import Element
from openmm import Vec3

from modules.docking_tasks import convert_to_pdbqt, dock
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

STANDARD_AA = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS",
    "ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"
}


# ---------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------

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


# ---------------------------------------------------------------------
# Core IFD system construction
# ---------------------------------------------------------------------

def build_ifd_system(protein_pdb: PDBFile, ligand_sdf: Path, config: dict):

    ff_path = config["induced_fit_docking"]["minimisation"]["protein_forcefield"]
    ligand_fixed = config["induced_fit_docking"]["minimisation"].get("ligand_fixed", True)

    # Start with protein
    modeller = Modeller(protein_pdb.topology, protein_pdb.positions)

    # Load ligand via OpenFF
    ligand_mol = Molecule.from_file(str(ligand_sdf))
    ligand_mol.name = "LIG"  # optional but nice

    if not ligand_mol.conformers:
        raise RuntimeError(f"Ligand has no 3D conformers: {ligand_sdf}")

    # Convert OpenFF → OpenMM topology
    ligand_top = ligand_mol.to_topology().to_openmm()

    # Ensure residue is named LIG
    for residue in ligand_top.residues():
        residue.name = "LIG"

    # GNINA SDF coordinates (numpy array, Å)
    coords = ligand_mol.conformers[0].to_openmm()
    ligand_positions = coords.in_units_of(unit.nanometer)

    # Add ligand ONCE
    modeller.add(ligand_top, ligand_positions)

    lig_atoms = [a for a in modeller.topology.atoms() if a.residue.name == "LIG"]
    lig_bonds = [
        b for b in modeller.topology.bonds()
        if b[0].residue.name == "LIG" or b[1].residue.name == "LIG"
    ]

    # print ligand info - debug

    logger.info(f"Ligand atoms: {len(lig_atoms)}, bonds: {len(lig_bonds)}")

    if len(lig_atoms) == 0 or len(lig_bonds) == 0:
        raise RuntimeError("Ligand was not correctly added to the topology")

    # Build system
    system_generator = SystemGenerator(
        forcefields=[ff_path],
        small_molecule_forcefield="openff-2.1.0",
        molecules=[ligand_mol],
    )

    system = system_generator.create_system(modeller.topology)

    # Optional: restrain ligand
    if ligand_fixed:
        ligand_atom_indices = [a.index for a in modeller.topology.atoms() if a.residue.name == "LIG"]
        add_positional_restraints(system, modeller, ligand_atom_indices, k=1000.0)

    return system, modeller


# ---------------------------------------------------------------------
# Thread-safe IFD task
# ---------------------------------------------------------------------

@register_task(
    "induced_fit_docking",
    category="Docking",
    description="Dock, minimise nearby residues and re-dock (thread-safe).",
)
def induced_fit_docking(backend, ligand, config, **kwargs):
    """
    Minimises the receptor around a docked ligand and then re-docks the ligand.
    Uses OpenFF SDF directly, no PDB conversion needed.
    """
    # ---------------------------
    # Per-ligand working directory
    # ---------------------------
    base_output = Path(config["output_dir"])
    ligand_dir = base_output / ligand["name"]
    ligand_dir.mkdir(parents=True, exist_ok=True)

    # ---------------------------
    # Load receptor
    # ---------------------------
    receptor_pdbqt = backend.cache.get("receptor_pdbqt")
    if receptor_pdbqt is None:
        raise RuntimeError("Receptor PDBQT not found in cache.")

    receptor_clean_pdb = Path(str(receptor_pdbqt).replace(".pdbqt", "_protonated.pdb"))
    protein_pdb = load_protein(receptor_clean_pdb)

    # ---------------------------
    # Ensure ligand docked pose (SDF)
    # ---------------------------
    if not ligand.get("pdbqt_paths"):
        convert_to_pdbqt(backend, ligand, config)

    docked_sdf = base_output / f"{ligand['name']}_conf0_docked.sdf"
    if not docked_sdf.exists():
        raise FileNotFoundError(f"Docked SDF not found: {docked_sdf}")

    # ---------------------------
    # Build & minimise system
    # ---------------------------
    system, modeller = build_ifd_system(protein_pdb, docked_sdf, config)
    
    minim_cfg = config["induced_fit_docking"]["minimisation"]
    tol = minim_cfg.get("tolerance", 0.5) # minim tolerance from the yaml

    integrator = openmm.LangevinIntegrator(
        300 * unit.kelvin,
        1 / unit.picosecond,
        0.002 * unit.picoseconds,
    )

    # Sanity checks

    assert system.getNumParticles() == modeller.topology.getNumAtoms()
    assert len(modeller.positions) == modeller.topology.getNumAtoms()

    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)

    simulation.minimizeEnergy(
        tolerance = tol * unit.kilojoule_per_mole / unit.nanometer,
        maxIterations=minim_cfg.get("max_steps", 500),
    )

    # ---------------------------
    # Write minimised receptor
    # ---------------------------
    minimised_pdb = ligand_dir / "protein_minimised.pdb"
    with open(minimised_pdb, "w") as f:
        PDBFile.writeFile(
            simulation.topology,
            simulation.context.getState(getPositions=True).getPositions(),
            f,
        )
    logger.info(f"[{ligand['name']}] Minimised receptor written")

    # ---------------------------
    # Re-dock using a LOCAL backend view
    # ---------------------------
    local_config = dict(config)
    local_config["protein"] = dict(config["protein"])
    local_config["protein"]["pdb_path"] = str(minimised_pdb)

    logger.info(f"[{ligand['name']}] Re-docking to minimised receptor")
    return dock(backend, ligand, local_config, **kwargs)

