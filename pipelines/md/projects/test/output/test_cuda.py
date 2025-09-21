from openmm.app import PDBFile

pdb = PDBFile('topology.pdb')
topology = pdb.topology

from openmm.app import ForceField

forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
system = forcefield.createSystem(topology)

from openmm import LangevinIntegrator
from openmm.unit import kelvin, picosecond, femtosecond

integrator = LangevinIntegrator(
    300*kelvin,
    1/picosecond,
    2*femtosecond
)

from openmm.app import Simulation
from openmm import Platform

platform = Platform.getPlatformByName('CUDA')
simulation = Simulation(topology, system, integrator, platform)

