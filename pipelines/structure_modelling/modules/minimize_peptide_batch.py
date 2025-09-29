# modules/minimize_peptides_openmm.py

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger
from pathlib import Path
import openmm.app as app
import openmm as mm
from openmm import unit

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_task("minimize_peptide_batch", category="Peptide modeling", description="Minimize peptide structures using OpenMM with optional implicit solvation")
def minimize_peptide_batch(backend, ligand, config, **kwargs):
    out_dir = Path(config["output_dir"])
    model_dict = backend.cache.get("peptide_models", {})
    minim_cfg = config.get("minimization", {})
    
    method = minim_cfg.get("method", "steepest_descent")
    steps = minim_cfg.get("steps", 5000)
    platform_name = minim_cfg.get("platform", "CPU")
    forcefield_name = minim_cfg.get("forcefield", "amber14-all.xml")
    implicit_solvent = minim_cfg.get("implicit_solvent", False)  # New flag
    
    logger.info(f"Minimizing {len(model_dict)} peptides with OpenMM on {platform_name} (implicit_solvent={implicit_solvent})...")

    platform = mm.Platform.getPlatformByName(platform_name)

    for name, pdb_path in model_dict.items():
        logger.info(f"Minimizing peptide: {name}")

        pdb = app.PDBFile(pdb_path)
        forcefield = app.ForceField(forcefield_name)
        
        modeller = app.Modeller(pdb.topology, pdb.positions)
        modeller.deleteWater()  # Optional: remove waters
        modeller.addHydrogens(forcefield)

        if implicit_solvent:
            system = forcefield.createSystem(
                modeller.topology,
                nonbondedMethod=app.NoCutoff,
                constraints=app.HBonds,
                implicitSolvent=app.OBC2)
        else:
            system = forcefield.createSystem(
                modeller.topology,
                nonbondedMethod=app.NoCutoff,
                constraints=app.HBonds)

        integrator = mm.LangevinIntegrator(
            300 * unit.kelvin,
            1 / unit.picosecond,
            0.002 * unit.picoseconds)
        
        simulation = app.Simulation(modeller.topology, system, integrator, platform)
        simulation.context.setPositions(modeller.positions)
        simulation.minimizeEnergy(maxIterations=steps)

        positions = simulation.context.getState(getPositions=True).getPositions()
        output_path = out_dir / name / "minimized.pdb"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            app.PDBFile.writeFile(simulation.topology, positions, f)
        
        logger.info(f"Minimized structure saved to: {output_path}")
