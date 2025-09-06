import os
import numpy as np
import yaml
from openmm.app import *
from openmm import *
from openmm.unit import *
from openmm import Platform
from pdbfixer import PDBFixer

class MDPipeline:
    def __init__(self, config: dict):
        self.config = config
        self.system = None
        self.simulation = None
        self.integrator = None
        self.topology = None
        self.positions = None
        self.restraints = []
        self.platform = Platform.getPlatformByName(config.get('platform', 'CPU'))

    @classmethod
    def from_yaml(cls, yaml_file: str):
        with open(yaml_file, 'r') as f:
            config = yaml.safe_load(f)
        return cls(config)

    def prepare_system(self):
        input_pdb = self.config['input_pdb']
        fixer = PDBFixer(filename=input_pdb)

        if self.config.get('add_missing_residues', True):
            fixer.findMissingResidues()
        if self.config.get('add_missing_atoms', True):
            fixer.findMissingAtoms()
        if self.config.get('add_hydrogens', True):
            fixer.addMissingHydrogens(pH=self.config.get('pH', 7.0))
        if self.config.get('add_solvent', True):
            fixer.addSolvent(
                boxSize=Vec3(*self.config['solvent_box']).value_in_unit(nanometers),
                ionicStrength=self.config.get('ionic_strength', 0.15)*molar
            )

        self.topology = fixer.topology
        self.positions = fixer.positions

        forcefield_files = self.config.get('forcefield', ['amber14-all.xml', 'amber14/tip3pfb.xml'])
        self.forcefield = ForceField(*forcefield_files)
        self.system = self.forcefield.createSystem(
            self.topology,
            nonbondedMethod=PME,
            nonbondedCutoff=1.0*nanometer,
            constraints=HBonds
        )

    def apply_restraints(self, restraint_level=1.0):
        """Apply position restraints to heavy atoms."""
        force = CustomExternalForce("0.5 * k * ((x - x0)^2 + (y - y0)^2 + (z - z0)^2)")
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")
        force.addPerParticleParameter("k")

        for atom in self.topology.atoms():
            if atom.element.name != 'hydrogen':
                index = atom.index
                pos = self.positions[index]
                k = restraint_level * kilocalories_per_mole / angstroms**2
                force.addParticle(index, [pos.x, pos.y, pos.z, k])
        self.system.addForce(force)
        self.restraints.append(force)

    def setup_simulation(self):
        self.integrator = LangevinIntegrator(
            self.config.get('temperature', 300)*kelvin,
            self.config.get('friction', 1.0)/picosecond,
            self.config.get('timestep', 2.0)*femtoseconds
        )
        self.simulation = Simulation(self.topology, self.system, self.integrator, self.platform)
        self.simulation.context.setPositions(self.positions)

        self.simulation.minimizeEnergy()
        self.simulation.context.setVelocitiesToTemperature(0*kelvin)

    def run_heating_equilibration(self):
        heating = self.config.get('heating', {})
        n_steps = heating.get('n_steps', 6)
        target_temp = heating.get('target_temp', 300) * kelvin
        equil_steps_per_temp = heating.get('steps_per_temp', 5000)
        restraint_schedule = heating.get('restraint_schedule', np.linspace(10.0, 0.0, n_steps))

        for i in range(n_steps):
            temp = (target_temp / n_steps) * (i + 1)
            k = restraint_schedule[i]
            print(f"[Heating] Step {i+1}/{n_steps} — Temp: {temp}, Restraint k: {k}")

            # Remove previous restraint force (if any)
            if self.restraints:
                self.system.removeForce(self.system.getNumForces() - 1)
                self.restraints.pop()
            self.apply_restraints(k)

            self.integrator.setTemperature(temp)
            self.simulation.context.reinitialize(preserveState=True)
            self.simulation.step(equil_steps_per_temp)

        # Final equilibration
        if self.restraints:
            self.system.removeForce(self.system.getNumForces() - 1)
            self.restraints = []
            self.simulation.context.reinitialize(preserveState=True)

        equil_steps = self.config.get('equilibration', {}).get('steps', 10000)
        print(f"[Equilibration] Final equilibration for {equil_steps} steps")
        self.simulation.step(equil_steps)

    def run_production(self):
        prod_config = self.config['production']
        length_ns = prod_config['length_ns']
        timestep_fs = self.config.get('timestep', 2.0)
        steps = int((length_ns * 1e6) / timestep_fs)

        self.simulation.reporters.append(DCDReporter(prod_config['output_dcd'], 1000))
        self.simulation.reporters.append(StateDataReporter(
            prod_config['output_log'],
            1000, step=True, temperature=True, energy=True, progress=True, remainingTime=True, speed=True,
            totalSteps=steps, separator='\t'
        ))

        print(f"[Production] Running {steps} steps ({length_ns} ns)")
        self.simulation.step(steps)

    def run(self):
        self.prepare_system()
        self.setup_simulation()
        self.run_heating_equilibration()
        self.run_production()