from openmm.app import *
from openmm import *
from openmm.unit import *
import os
from openmm.unit import picoseconds

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

    @register_task("heat_and_equilibrate", category='Molecular dynamics', description="Heating and equilibration.")
    def heat_and_equilibrate(self):
        cfg = self.config["heat_and_equilibrate"]
        num_steps = cfg.get("num_heating_steps", 6)
        heating_increment = cfg.get("heating_increment", 298.15 / num_steps)
        target_temp = cfg.get("target_temp", 300)
        steps_per_round = cfg.get("steps_per_round", 5000)
        restraints = cfg.get("restraint_strengths", [5.0, 4.0, 3.0, 2.0, 1.0, 0.0])

        # Backbone atoms heating/equil constraints
        backbone_atoms = [atom.index for atom in self.topology.atoms() if atom.name in {"N", "CA", "C", "O"}]

        logger.info(f"Setting up {num_steps} heating/equilibration steps.")

        for i in range(num_steps):
            temp = heating_increment * (i + 1)
            restraint_k = restraints[i] if i < len(restraints) else 0

            logger.info(f"Heating step {i + 1}/{num_steps}: Temperature = {temp:.2f} K, Restraint strength = {restraint_k}")

            self.integrator.setTemperature(temp * kelvin)

            for idx in reversed(range(self.system.getNumForces())):
                force = self.system.getForce(idx)
                if isinstance(force, CustomExternalForce):
                    self.system.removeForce(idx)

            if restraint_k > 0:
                force = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
                force.addPerParticleParameter("x0")
                force.addPerParticleParameter("y0")
                force.addPerParticleParameter("z0")
                force.addGlobalParameter("k", restraint_k)

                state = self.simulation.context.getState(getPositions=True)
                positions = state.getPositions()

                for atom_idx in backbone_atoms:
                    pos = positions[atom_idx]
                    force.addParticle(atom_idx, [pos.x, pos.y, pos.z])

                self.system.addForce(force)
                positions = self.simulation.context.getState(getPositions=True).getPositions()
                self.simulation.context.reinitialize()
                self.simulation.context.setPositions(positions)

            self.simulation.step(steps_per_round)

    @register_task("production", category='Molecular dynamics', description="Run production simulation.")
    def production(self):
        cfg = self.config["production"]
        ns = cfg.get("length_ns", 1)
        timestep = self.integrator.getStepSize() 
        
        # Timestep -> ns
        timestep_ns = timestep.value_in_unit(picoseconds) / 1000.0

        steps_per_ns = int(1.0 / timestep_ns)
        total_steps = int(ns * steps_per_ns)

        output_trajectory = cfg.get("output_trajectory", "output")
        output_logfile = cfg.get("output_logfile", "poo.log")

        output_dir = os.path.dirname(output_trajectory)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        step_chunk = 1000  # Steps per update

        output_dir = os.path.dirname(output_trajectory)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Save checkpoint as "restart.chk" in the output directory
        checkpoint_file = os.path.join(output_dir, "restart.chk")

        self.simulation.reporters.append(DCDReporter(output_trajectory, step_chunk))
        self.simulation.reporters.append(StateDataReporter(output_logfile,
                                                        1000, step=True, temperature=True,
                                                        potentialEnergy=True,totalEnergy=True,
                                                        density=True,
                                                        progress=True, remainingTime=True,
                                                        speed=True, totalSteps=total_steps,
                                                        separator='\t'))

        logger.info(f"Running Production stage for {ns} ns")

        
        with tqdm(total=total_steps, desc="Production", unit="steps") as pbar:
            for _ in range(0, total_steps, step_chunk):
                steps_to_run = min(step_chunk, total_steps - self.simulation.currentStep)
                self.simulation.step(steps_to_run)

                # Save restart every 1000 steps.
                self.simulation.saveCheckpoint(checkpoint_file)

                logger.debug(f"Saving checkpoint to {checkpoint_file} every {steps_to_run} steps. Current step: {self.simulation.currentStep}")

                pbar.update(steps_to_run)

