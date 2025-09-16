from openff.toolkit.topology import Molecule
from simtk.openmm.app import PDBFile
from simtk.openmm import Vec3
from simtk.unit import angstrom
from openmm.app.element import get_by_symbol
from openmm.app.topology import Topology
from openmm.app import Modeller

def load_ligand_from_sdf(sdf_path):
    mol = Molecule.from_file(sdf_path, file_format='sdf')
    mol.generate_conformers(n_conformers=1)
    pdb_path = sdf_path.replace(".sdf", "_ligand.pdb")
    mol.to_file(pdb_path, "pdb")
    return pdb_path
