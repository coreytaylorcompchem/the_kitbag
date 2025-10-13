import numpy as np
import torch
from torch_geometric.data import Data, Batch
from rdkit import Chem
from rdkit.Chem import Descriptors, rdPartialCharges, rdMolDescriptors

# ==========================================================
# Custom Data structure
# ==========================================================
class MoleculeData(Data):
    """Custom PyG Data object to prevent global features from being concatenated improperly."""
    def __inc__(self, key, value, *args, **kwargs):
        return super().__inc__(key, value)

    def __cat_dim__(self, key, value, *args, **kwargs):
        if key == 'global_features':
            return None  # stack, don’t concatenate
        return super().__cat_dim__(key, value)

# ==========================================================
# Helper functions
# ==========================================================
def safe_value(val):
    try:
        if np.isnan(val) or np.isinf(val):
            return 0.0
        return float(val)
    except Exception:
        return 0.0

def one_hot_encoding(x, allowable_set):
    if x not in allowable_set:
        x = allowable_set[-1]
    return [int(x == s) for s in allowable_set]

# ==========================================================
# Periodic table info
# ==========================================================
electronegativity_dict = {
    1: 2.20, 6: 2.55, 7: 3.04, 8: 3.44, 9: 3.98, 15: 2.19,
    16: 2.58, 17: 3.16, 35: 2.96, 53: 2.66,
}

metals = {
    3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27,
    28, 29, 30, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
    48, 49, 55, 56, 57, 72, 73, 74, 75, 76, 77, 78, 79, 80,
    81, 82, 83, 87, 88, 89, 104, 105, 106, 107, 108, 109,
    110, 111, 112, 113, 114, 115, 116
}

# ==========================================================
# Atom & Bond features
# ==========================================================
def get_atom_features(atom, mol=None):
    pt = Chem.GetPeriodicTable()
    features = []

    features += one_hot_encoding(atom.GetSymbol(),
                                ['C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I', 'H',
                                 'B', 'Si', 'Se', 'As', 'Al', 'Zn', 'Cu', 'Ni', 'Fe', 'other'])
    features += one_hot_encoding(atom.GetHybridization(),
                                [Chem.rdchem.HybridizationType.SP,
                                 Chem.rdchem.HybridizationType.SP2,
                                 Chem.rdchem.HybridizationType.SP3,
                                 Chem.rdchem.HybridizationType.SP3D,
                                 Chem.rdchem.HybridizationType.SP3D2,
                                 Chem.rdchem.HybridizationType.UNSPECIFIED])

    features += [
        atom.GetDegree(),
        atom.GetFormalCharge(),
        atom.GetNumRadicalElectrons(),
        int(atom.GetIsAromatic())
    ]

    atomic_num = atom.GetAtomicNum()
    features += [
        pt.GetAtomicWeight(atomic_num) / 200,
        pt.GetRvdw(atomic_num) / 2.5,
        electronegativity_dict.get(atomic_num, 0.0) / 4.0,
        pt.GetNOuterElecs(atomic_num) / 8.0,
        1.0 if atomic_num in metals else 0.0
    ]

    if mol is not None:
        try:
            charge = float(atom.GetProp('_GasteigerCharge'))
            charge = safe_value(charge)
        except Exception:
            charge = 0.0
        features += [charge]
    else:
        features += [0.0]

    return np.array(features, dtype=np.float32)

def get_bond_features(bond):
    bond_type = bond.GetBondType()
    return np.array([
        int(bond_type == Chem.rdchem.BondType.SINGLE),
        int(bond_type == Chem.rdchem.BondType.DOUBLE),
        int(bond_type == Chem.rdchem.BondType.TRIPLE),
        int(bond_type == Chem.rdchem.BondType.AROMATIC),
        int(bond.GetIsConjugated()),
        int(bond.IsInRing())
    ], dtype=np.float32)

# ==========================================================
# Global molecular descriptors
# ==========================================================
descriptor_functions = [
    Descriptors.MolWt, Descriptors.MolLogP, Descriptors.NumHDonors,
    Descriptors.NumHAcceptors, Descriptors.TPSA, Descriptors.FractionCSP3,
    Descriptors.HeavyAtomCount, Descriptors.NumRotatableBonds, Descriptors.RingCount,
    Descriptors.NumValenceElectrons, Descriptors.NumHeteroatoms,
    rdMolDescriptors.CalcNumAliphaticRings, rdMolDescriptors.CalcNumAromaticRings,
    rdMolDescriptors.CalcNumSaturatedRings, rdMolDescriptors.CalcLabuteASA,
    lambda mol: rdMolDescriptors.CalcCrippenDescriptors(mol)[0],
    lambda mol: rdMolDescriptors.CalcCrippenDescriptors(mol)[1],
]

# ==========================================================
# Graph conversion
# ==========================================================
def mol_to_graph(smiles: str, label: float = None):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    try:
        rdPartialCharges.ComputeGasteigerCharges(mol)
    except Exception:
        pass

    x = torch.tensor([get_atom_features(a, mol) for a in mol.GetAtoms()], dtype=torch.float)
    edge_index, edge_attr = [], []

    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bfeat = get_bond_features(bond)
        edge_index += [[i, j], [j, i]]
        edge_attr += [bfeat, bfeat]

    edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    edge_attr = torch.tensor(edge_attr, dtype=torch.float)

    global_features = []
    for func in descriptor_functions:
        try:
            global_features.append(safe_value(func(mol)))
        except Exception:
            global_features.append(0.0)
    global_features = torch.tensor(global_features, dtype=torch.float32)

    data = MoleculeData(
        x=x,
        edge_index=edge_index,
        edge_attr=edge_attr,
        global_features=global_features
    )
    if label is not None:
        data.y = torch.tensor([label], dtype=torch.float)
    data.smiles = smiles
    return data

def compute_global_features(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return torch.zeros(len(descriptor_functions), dtype=torch.float)
    feats = [safe_value(fn(mol)) for fn in descriptor_functions]
    return torch.tensor(feats, dtype=torch.float)

def custom_collate(batch_list):
    return Batch.from_data_list(batch_list)
