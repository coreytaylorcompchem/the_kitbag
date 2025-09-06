from md import MDPipeline
from openmmml import MLPotential
from openmm import *

class QMMMPipeline(MDPipeline):
    def setup_simulation(self):
        # Add ML/MM force to system (ANI2x for QM region)
        print("[QMMM] Using ANI ML potential for QMMM region")
        mlpotential = MLPotential("ANI2x")
        qm_atoms = self.config['qmmm']['qm_atoms']  # list of atom indices
        ml_system = mlpotential.createPotential(self.topology, qm_atoms)
        self.system.addForce(ml_system)
        
        # Now call normal setup
        super().setup_simulation()
