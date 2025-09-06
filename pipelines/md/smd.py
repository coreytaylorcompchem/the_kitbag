from md import MDPipeline
from openmm import *
from openmm.unit import *

class SMDPipeline(MDPipeline):
    def add_smd_force(self):
        smd_config = self.config['smd']
        atom1 = smd_config['atom1']  # index
        atom2 = smd_config['atom2']  # index
        k = smd_config.get('spring_constant', 5.0) * kilocalories_per_mole / angstroms**2
        velocity = smd_config.get('velocity', 0.001) * angstroms / picoseconds

        # Create custom bond force
        force = CustomBondForce("0.5 * k * (r - r0(t))^2")
        force.addPerBondParameter("k")

        # Time-dependent r0
        r0 = smd_config.get("initial_distance", 1.0) * nanometers
        expression = f"{r0.value_in_unit(nanometers)} + {velocity.value_in_unit(nanometers/picoseconds)} * t"
        force.setGlobalParameterDefaultValue("t", 0)
        force.addGlobalParameter("t", 0)
        force.addBond(atom1, atom2, [k])
        self.system.addForce(force)

        # Store the force for updating r0 with time
        self.smd_force = force

    def run_production(self):
        self.add_smd_force()

        prod_config = self.config['production']
        steps = int((prod_config['length_ns'] * 1e6) / self.config['timestep'])

        # Add reporters
        self.simulation.reporters.append(DCDReporter(prod_config['output_dcd'], 1000))
        self.simulation.reporters.append(StateDataReporter(prod_config['output_log'], 1000,
                                                           step=True, potentialEnergy=True, temperature=True))

        # Run with time update
        print(f"[SMD] Running {steps} steps")
        for step in range(steps):
            time_ps = step * self.config['timestep'] * femtoseconds
            self.simulation.context.setParameter("t", time_ps.value_in_unit(picoseconds))
            self.simulation.step(1)
