import os

from openmm.app import *
from openmm import *
from openmm.unit import *
from openmm.unit import picoseconds
from openmm import Vec3
from openmm.unit import nanometer, kelvin, atmosphere, picoseconds

import numpy as np
from tqdm import tqdm

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

class MDWorkflow:
    def __init__(self, backend, config):
        self.backend = backend
        self.config = config
        self.integrator = None
        self.simulation = None

    @register_task("prepare_system", category='Molecular dynamics',
                   description="Load inputs, cap chains, parameterise ligand and protein, solvate and save final topology.")
    def prepare_system(self):
        self.backend.prepare_system(self.config["prepare_system"])
        self.system = self.backend.system
        self.topology = self.backend.topology
        self.positions = self.backend.positions

    @register_task("setup_simulation", category='Molecular dynamics', description="Set up integrator and simulation.")
    def setup_simulation(self):
        cfg = self.config["setup_simulation"]
        self.integrator = LangevinIntegrator(
            cfg["temperature"]*kelvin,
            1.0/picoseconds,
            cfg["timestep"]*femtoseconds
        )
        self.simulation = Simulation(self.topology, self.system, self.integrator,
                                     Platform.getPlatformByName(cfg["platform"]))
        self.simulation.context.setPositions(self.positions)

    @register_task("minimize", category='Molecular dynamics', description="Initial energy minimization.")
    def minimize(self):
        self.simulation.minimizeEnergy()

    @register_task("heat_and_equilibrate", category='Molecular dynamics',
            description="Heating and equilibration.")
    def heat_and_equilibrate(self):
        cfg = self.config["heat_and_equilibrate"]

        # -------------------------------
        # 1. HEATING PHASE
        # -------------------------------
        heating_cfg = cfg.get("heating", {})
        num_steps = heating_cfg.get("num_steps", 6)
        steps_per_round = heating_cfg.get("steps_per_round", 5000)
        target_temp = heating_cfg.get("target_temp", 300)
        restraints = heating_cfg.get("restraint_strengths", [10.0, 5.0, 2.0, 1.0, 0.5, 0.0])

        backbone_atoms = [atom.index for atom in self.topology.atoms() if atom.name in {"N","CA","C","O"}]
        heating_increment = target_temp / num_steps

        logger.info(f"Starting heating: {num_steps} rounds → {target_temp} K")

        # Need to create a single CustomExternalForce for positional restraints
        force = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")
        force.addGlobalParameter("k", restraints[0])  # initial max restraint

        # Get initial positions
        state = self.simulation.context.getState(getPositions=True)
        positions = state.getPositions()

        for atom_idx in backbone_atoms:
            pos = positions[atom_idx]
            force.addParticle(atom_idx, [pos.x, pos.y, pos.z])

        self.system.addForce(force)
        self.simulation.context.reinitialize(preserveState=True)
        self.simulation.context.setPositions(positions)

        # Heating loop
        for i in range(num_steps):
            temp = heating_increment * (i + 1)
            k_value = restraints[i] if i < len(restraints) else 0.0
            logger.info(f"Heating round {i+1}/{num_steps}: T = {temp:.2f} K, restraint = {k_value}")

            # Update temperature and restraint
            self.integrator.setTemperature(temp * kelvin)
            self.simulation.context.setParameter("k", k_value)

            self.simulation.step(steps_per_round)

        # Remove restraints after heating
        self.system.removeForce(self.system.getNumForces() - 1)
        self.simulation.context.reinitialize(preserveState=True)
        self.simulation.context.setPositions(self.simulation.context.getState(getPositions=True).getPositions())

        # -------------------------------
        # 2. EQUILIBRATION PHASE
        # -------------------------------
        equil_cfg = cfg.get("equilibration", {})
        ensemble = equil_cfg.get("ensemble", "NPT").upper()
        steps = equil_cfg.get("steps", 100000)
        equil_temp = equil_cfg.get("target_temp", target_temp)

        logger.info(f"Starting equilibration: {ensemble} ensemble for {steps} steps at {equil_temp} K")

        # Add barostat if NPT and not already present
        has_barostat = any(isinstance(f, MonteCarloBarostat) for f in
                        [self.system.getForce(i) for i in range(self.system.getNumForces())])
        if ensemble == "NPT" and not has_barostat:
            pressure = equil_cfg.get("pressure", 1.0) * atmosphere
            barostat_cfg = equil_cfg.get("barostat", {})
            frequency = barostat_cfg.get("frequency", 25)
            barostat = MonteCarloBarostat(pressure, equil_temp * kelvin, frequency)
            self.system.addForce(barostat)
            self.simulation.context.reinitialize(preserveState=True)
            self.simulation.context.setPositions(positions)
            logger.info(f"Added MonteCarloBarostat: {pressure} atm, {equil_temp} K")

        # Ensure integrator temperature
        self.integrator.setTemperature(equil_temp * kelvin)

        # Run equilibration
        step_chunk = 1000
        with tqdm(total=steps, desc="Equilibration", unit="steps") as pbar:
            current_step = 0
            while current_step < steps:
                n = min(step_chunk, steps - current_step)
                self.simulation.step(n)
                current_step += n
                pbar.update(n)

        # Diagnostics
        state = self.simulation.context.getState(getEnergy=True, getPositions=True)
        box_vectors = state.getPeriodicBoxVectors(asNumpy=True)
        volume_nm3 = (box_vectors[0][0] * box_vectors[1][1] * box_vectors[2][2]).value_in_unit(nanometer**3)

        mass = sum([self.system.getParticleMass(i) for i in range(self.system.getNumParticles())])
        mass_g = mass.value_in_unit(dalton) * 1.66054e-24
        volume_mL = volume_nm3 * 1e-21
        density = mass_g / volume_mL

        logger.info(f"Equilibration complete. Box volume: {volume_nm3:.3f} nm³, density: {density:.3f} g/mL")
        if abs(density - 1.0) > 0.1:
            logger.warning("Density deviates from expected ~1.0 g/mL. Check for bubbles or equilibration issues.")


    @register_task("production", category='Molecular dynamics', description="Run production simulation.")
    def production(self):

        cfg = self.config["production"]
        ns = cfg.get("length_ns", 1)
        ensemble = cfg.get("ensemble", "NPT").upper()
        target_temp = cfg.get("target_temp", self.integrator.getTemperature().value_in_unit(kelvin))
        pressure = cfg.get("pressure", 1.0)  # atm

        # Get current positions, velocities, and box vectors
        state = self.simulation.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
        positions = state.getPositions()
        velocities = state.getVelocities()
        box_vectors = state.getPeriodicBoxVectors()

        #### To center in the box 
        # Convert positions to plain floats in nm
        coords = np.array([
            [
                p.x.value_in_unit(nanometer) if hasattr(p.x, "value_in_unit") else float(p.x),
                p.y.value_in_unit(nanometer) if hasattr(p.y, "value_in_unit") else float(p.y),
                p.z.value_in_unit(nanometer) if hasattr(p.z, "value_in_unit") else float(p.z)
            ]
            for p in positions
        ], dtype=float)

        # Get backbone atom indices
        backbone_atoms = [atom.index for atom in self.topology.atoms() if atom.name in {"N", "CA", "C", "O"}]

        # Compute Center of Mass
        com = np.mean(coords[backbone_atoms], axis=0)

        # Box lengths
        box_lengths = np.array([
            box_vectors[0][0].value_in_unit(nanometer) if hasattr(box_vectors[0][0], "value_in_unit") else float(box_vectors[0][0]),
            box_vectors[1][1].value_in_unit(nanometer) if hasattr(box_vectors[1][1], "value_in_unit") else float(box_vectors[1][1]),
            box_vectors[2][2].value_in_unit(nanometer) if hasattr(box_vectors[2][2], "value_in_unit") else float(box_vectors[2][2])
        ], dtype=float)

        # Compute shift
        shift = box_lengths / 2 - com

        # Apply shift to positions
        for i in range(len(positions)):
            x, y, z = coords[i] + shift
            positions[i] = Vec3(x, y, z) * nanometer

        # Update simulation context
        self.simulation.context.setPositions(positions)
        self.simulation.context.setVelocities(velocities)
        self.simulation.context.setPeriodicBoxVectors(*box_vectors)
        self.integrator.setTemperature(target_temp * kelvin)
        self.simulation.currentStep = 0

        # Attach barostat
        has_barostat = any(isinstance(f, MonteCarloBarostat) for f in
                        [self.system.getForce(i) for i in range(self.system.getNumForces())])
        if ensemble == "NPT" and not has_barostat:
            barostat = MonteCarloBarostat(pressure * atmosphere, target_temp * kelvin, 25)
            self.system.addForce(barostat)
            self.simulation.context.reinitialize(preserveState=True)

            # Restore positions and velocities
            self.simulation.context.setPositions(positions)
            self.simulation.context.setVelocities(velocities)
            self.simulation.context.setPeriodicBoxVectors(*box_vectors)
            logger.info(f"Added MonteCarloBarostat: {pressure} atm, {target_temp} K")

        # Compute total steps
        timestep = self.integrator.getStepSize()
        timestep_ns = timestep.value_in_unit(picoseconds) / 1000.0
        steps_per_ns = int(1.0 / timestep_ns)
        total_steps = int(ns * steps_per_ns)

        # Setup outputs
        output_trajectory = cfg.get("output_trajectory", "output.dcd")
        output_logfile = cfg.get("output_logfile", "output.log")
        output_dir = os.path.dirname(output_trajectory)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        step_chunk = 1000 #TODO add this to yaml
        checkpoint_file = os.path.join(output_dir, "restart.chk")

        self.simulation.reporters.append(DCDReporter(output_trajectory, step_chunk))
        self.simulation.reporters.append(StateDataReporter(
            output_logfile, step_chunk,
            step=True, temperature=True, potentialEnergy=True, totalEnergy=True,
            density=True, progress=True, remainingTime=True, speed=True,
            totalSteps=total_steps, separator='\t'
        ))

        # Run prod
        logger.info(f"Running Production: {ns} ns ({total_steps} steps) in {ensemble} ensemble.")

        with tqdm(total=total_steps, desc="Production", unit="steps") as pbar:
            while self.simulation.currentStep < total_steps:
                steps_to_run = min(step_chunk, total_steps - self.simulation.currentStep)
                self.simulation.step(steps_to_run)
                self.simulation.saveCheckpoint(checkpoint_file)
                pbar.update(steps_to_run)

                # # Wrap atoms every 10k steps
                # This doesn't work and may not be necessary anyway so commenting out for now.
                # if self.simulation.currentStep % 10000 == 0:
                #     state = self.simulation.context.getState(getPositions=True)
                #     positions = state.getPositions()
                #     for i in range(len(positions)):
                #         x = (positions[i].x / nanometer) % box_lengths[0]
                #         y = (positions[i].y / nanometer) % box_lengths[1]
                #         z = (positions[i].z / nanometer) % box_lengths[2]
                #         positions[i] = Vec3(float(x), float(y), float(z)) * nanometer
                #     self.simulation.context.setPositions(positions)





