from openmm.app import *
from openmm import *
from openmm.unit import *
import os

class MDWorkflow:
    def __init__(self, backend, config):
        self.backend = backend
        self.config = config
        self.integrator = None
        self.simulation = None

    def prepare_system(self):
        # Call the backend method that actually prepares the system
        self.backend.prepare_system()
        self.system = self.backend.system
        self.topology = self.backend.topology
        self.positions = self.backend.positions

    def setup_simulation(self):
        self.integrator = LangevinIntegrator(
            self.config["simulation"]["temperature"]*kelvin,
            1.0/picosecond,
            self.config["simulation"]["timestep"]*femtoseconds
        )
        self.simulation = Simulation(self.topology, self.system, self.integrator, self.backend.platform)
        self.simulation.context.setPositions(self.positions)

    def minimize(self):
        self.simulation.minimizeEnergy()

    def heat_and_equilibrate(self):
        num_steps = self.config["equilibration"]["num_heating_steps"]
        increment = self.config["equilibration"]["heating_increment"]
        restraints = self.config["equilibration"]["restraint_strengths"]

        for i in range(num_steps):
            temp = increment * (i + 1)
            k = restraints[i] if i < len(restraints) else 0

            self.integrator.setTemperature(temp * kelvin)

            if k > 0:
                force = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
                force.addPerParticleParameter("x0")
                force.addPerParticleParameter("y0")
                force.addPerParticleParameter("z0")
                force.addGlobalParameter("k", k)
                for idx, atom in enumerate(self.topology.atoms()):
                    pos = self.simulation.context.getState(getPositions=True).getPositions()[idx]
                    force.addParticle(idx, [pos.x, pos.y, pos.z])
                self.system.addForce(force)

            self.simulation.step(5000)

    def run_production(self):
        ns = self.config["simulation"]["total_ns"]
        steps = int(ns * self.config["simulation"]["steps_per_ns"])

        self.simulation.reporters.append(DCDReporter(self.config["output"]["trajectory"], 1000))
        self.simulation.reporters.append(StateDataReporter(self.config["output"]["logfile"],
                                                           1000, step=True, temperature=True,
                                                           progress=True, remainingTime=True,
                                                           speed=True, totalSteps=steps,
                                                           separator='\t'))
        self.simulation.step(steps)
