from backends.xtb import XTBBackend
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdMolTransforms
from sklearn.decomposition import PCA
import numpy as np

def get_backend_from_config(task_name, config):
    """
    Given a task name and the global config, return the appropriate backend instance.
    """
    task_cfg = config.get(task_name, {})
    backend_name = task_cfg.get("backend")

    if backend_name == "xtb":
        return XTBBackend()
    else:
        raise ValueError(f"Unknown backend '{backend_name}' for task '{task_name}'")

def calc_rmsd(mol1, mol2):
    conf1 = mol1.GetConformer()
    conf2 = mol2.GetConformer()
    pos1 = np.array(conf1.GetPositions())
    pos2 = np.array(conf2.GetPositions())
    diff = pos1 - pos2
    return np.sqrt((diff ** 2).mean())

def is_rotor_bond(bond):
    if bond.IsInRing():
        return False
    if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
        return False
    begin_atom = bond.GetBeginAtom()
    end_atom = bond.GetEndAtom()
    # Avoid terminal bonds (e.g., bonds to hydrogens or end groups)
    if begin_atom.GetDegree() <= 1 or end_atom.GetDegree() <= 1:
        return False
    return True
