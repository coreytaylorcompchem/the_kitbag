from md import MDPipeline
from openmm.app import Metadynamics
from openmm import *

class MetaDynamicsPipeline(MDPipeline):
    def setup_simulation(self):
        super().setup_simulation()

        # Define CV — distance between two atoms
        cv_config = self.config['metadynamics']
        atom1 = cv_config['atom1']
        atom2 = cv_config['atom2']
        self.collective_variable = CustomBondForce("r")
        self.collective_variable.addBond(atom1, atom2)
        self.cv_index = self.system.addForce(self.collective_variable)

    def run_production(self):
        cv = self.collective_variable

        # Setup MetaD
        bias = Metadynamics(
            self.system, [cv],
            temperature=self.config['temperature'] * kelvin,
            biasFactor=10.0,
            height=1.0 * kilojoule_per_mole,
            frequency=500,
            saveFrequency=1000,
            biasDir='metad_output'
        )

        prod_config = self.config['production']
        steps = int((prod_config['length_ns'] * 1e6) / self.config['timestep'])

        self.simulation.reporters.append(DCDReporter(prod_config['output_dcd'], 1000))
        self.simulation.reporters.append(StateDataReporter(prod_config['output_log'], 1000, temperature=True))

        print(f"[MetaD] Running {steps} steps with bias")
        bias.step(self.simulation, steps)
