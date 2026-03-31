import numpy as np

from rdkit.Chem import AllChem
from rdkit import DataStructs

def compute_morgan_fp(mol, n_bits=512):
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=n_bits)
    arr = np.zeros((n_bits,), dtype=np.float32)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr