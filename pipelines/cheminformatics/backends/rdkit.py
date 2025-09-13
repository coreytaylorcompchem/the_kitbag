from rdkit import Chem
from backends.base_backend import BaseBackend

class RDKitBackend(BaseBackend):
    def load_molecules(self, source):
        suppl = Chem.SDMolSupplier(source)
        return [mol for mol in suppl if mol is not None]

    def compute_descriptors(self, mol):
        return {
            'MolWt': Chem.Descriptors.MolWt(mol),
            'LogP': Chem.Descriptors.MolLogP(mol),
        }

    def filter_molecule(self, mol, mw_cutoff=500):
        return Chem.Descriptors.MolWt(mol) < mw_cutoff
