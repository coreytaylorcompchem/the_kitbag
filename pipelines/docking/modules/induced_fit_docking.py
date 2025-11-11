import os
from pathlib import Path
import tempfile
import shutil

import openmm
import openmm.app as app
import openmm.unit as unit
from openmm.app import PDBFile, Modeller, ForceField, Simulation
from openmm.app import PDBReporter

from modules.docking_tasks import convert_to_pdbqt, dock

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def load_protein(pdb_path: Path):
    return PDBFile(str(pdb_path))

def apply_induced_fit_modifications(protein_pdb: PDBFile, ligand_pdb: Path, config: dict, ligand_fixed: bool):
    """
    Builds an OpenMM system for minimisation of selected protein residues.
    """
    ff_path = config["induced_fit_docking"]["minimisation"]["protein_forcefield"]
    forcefield = ForceField(ff_path)
    
    modeller = Modeller(protein_pdb.topology, protein_pdb.positions)
    
    # Add ligand to modeller
    if ligand_pdb.exists():
        ligand = PDBFile(str(ligand_pdb))
        modeller.add(ligand.topology, ligand.positions)
    else:
        logger.warning(f"Ligand PDB not found at {ligand_pdb}, continuing with protein only.")

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=None
    )

    # If ligand_fixed, zero out forces for ligand atoms
    if ligand_fixed and ligand_pdb.exists():
        ligand_atoms = [atom.index for atom in modeller.topology.atoms() if atom.residue.name not in ("ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL")]
        for idx in ligand_atoms:
            for f in system.getForces():
                try:
                    f.setParticleMass(idx, 0.0*unit.dalton)
                except Exception:
                    pass

    return system, modeller

@register_task("induced_fit_docking",
               category="Docking",
               description="Dock, minimise nearby residues and re-dock.")
def induced_fit_docking(backend, ligand, config, **kwargs):
    """
    Minimises the receptor around a docked ligand and then re-docks the ligand.
    """
    output_dir = Path(config.get("output_dir", "output"))
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load receptor and ligand
    receptor_pdbqt = backend.cache.get("receptor_pdbqt")
    if receptor_pdbqt is None:
        raise RuntimeError("Receptor PDBQT not found in cache. Run 'prepare_receptor_pdbqt' first.")
    
    receptor_clean_pdb = Path(str(receptor_pdbqt).replace(".pdbqt", "_clean.pdb"))
    ligand_pdbqt_paths = ligand.get("pdbqt_paths", [])
    
    if not ligand_pdbqt_paths:
        logger.warning(f"No PDBQT paths found for ligand {ligand['name']}, generating default conformer PDBQT first.")
        convert_to_pdbqt(backend, ligand, config)
        ligand_pdbqt_paths = ligand.get("pdbqt_paths", [])

    # --- Step 1: Minimise protein around ligand ---
    ligand_fixed = config["induced_fit_docking"]["minimisation"].get("ligand_fixed", True)
    protein_pdb = load_protein(receptor_clean_pdb)
    
    # Use first conformer for minimisation
    ligand_pdb = Path(str(ligand_pdbqt_paths[0]).replace(".pdbqt", ".pdb"))
    system, modeller = apply_induced_fit_modifications(protein_pdb, ligand_pdb, config, ligand_fixed)

    integrator = openmm.LangevinIntegrator(
        300*unit.kelvin,
        1/unit.picosecond,
        0.002*unit.picoseconds
    )

    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)

    max_steps = config["induced_fit_docking"]["minimisation"].get("max_steps", 500)
    tol = config["induced_fit_docking"]["minimisation"].get("tolerance", 1e-4)

    logger.info(f"Starting protein minimisation: max_steps={max_steps}, tolerance={tol}")
    simulation.minimizeEnergy(tolerance=tol*unit.kilojoule/unit.mole, maxIterations=max_steps)
    logger.info("Protein minimisation complete.")

    # Save minimised receptor
    minimised_pdb = output_dir / f"{ligand['name']}_protein_minimised.pdb"
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), open(minimised_pdb, 'w'))
    logger.info(f"Minimised receptor saved at {minimised_pdb}")
    backend.cache["receptor_pdbqt"] = minimised_pdb  # overwrite cached receptor for next docking

    # --- Step 2: Re-dock ligand to minimised protein ---
    logger.info(f"Re-docking ligand {ligand['name']} to minimised protein...")
    docking_outputs = dock(backend, ligand, config, **kwargs)

    return docking_outputs
