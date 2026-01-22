import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem import Descriptors, rdPartialCharges
from torch_geometric.data import Data

# ---------- Data class ----------

class MoleculeData(Data):
    def __cat_dim__(self, key, value, *args, **kwargs):
        if key == 'global_features':
            return None
        return super().__cat_dim__(key, value)

# ---------- Descriptor definitions ----------

descriptor_functions = [
    Descriptors.MolWt,
    Descriptors.MolLogP,
    Descriptors.NumHDonors,
    Descriptors.NumHAcceptors,
    Descriptors.TPSA,
    Descriptors.FractionCSP3,
    Descriptors.HeavyAtomCount,
    Descriptors.NumRotatableBonds,
    Descriptors.RingCount
]

electronegativity_dict = {
    1: 2.20, 6: 2.55, 7: 3.04, 8: 3.44, 9: 3.98,
    15: 2.19, 16: 2.58, 17: 3.16, 35: 2.96, 53: 2.66,
}

metals = {
    3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 30, 37, 38, 39, 40, 41, 42,
    43, 44, 45, 46, 47, 48, 49, 55, 56, 57, 72,
    73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83,
    87, 88, 89
}

# ---------- Utility ----------

def safe_value(val):
    try:
        if np.isnan(val) or np.isinf(val):
            return 0.0
        return float(val)
    except:
        return 0.0

def one_hot_encoding(x, allowable):
    if x not in allowable:
        x = allowable[-1]
    return [int(x == s) for s in allowable]

# ---------- Atom and bond features ----------

def get_atom_features(atom, mol=None):
    pt = Chem.GetPeriodicTable()
    feats = []

    feats += one_hot_encoding(atom.GetSymbol(),
        ['C','N','O','F','P','S','Cl','Br','I','H',
         'B','Si','Se','As','Al','Zn','Cu','Ni','Fe','other'])

    feats += one_hot_encoding(atom.GetHybridization(), [
        Chem.rdchem.HybridizationType.SP,
        Chem.rdchem.HybridizationType.SP2,
        Chem.rdchem.HybridizationType.SP3,
        Chem.rdchem.HybridizationType.SP3D,
        Chem.rdchem.HybridizationType.SP3D2,
        Chem.rdchem.HybridizationType.UNSPECIFIED
    ])

    feats += [
        atom.GetDegree(),
        atom.GetFormalCharge(),
        atom.GetNumRadicalElectrons(),
        int(atom.GetIsAromatic())
    ]

    Z = atom.GetAtomicNum()
    feats += [
        pt.GetAtomicWeight(Z) / 200,
        pt.GetRvdw(Z) / 2.5,
        electronegativity_dict.get(Z, 0.0) / 4.0,
        pt.GetNOuterElecs(Z) / 8.0,
        1.0 if Z in metals else 0.0
    ]

    if mol is not None:
        try:
            charge = safe_value(atom.GetProp('_GasteigerCharge'))
        except:
            charge = 0.0
        feats += [charge]
    else:
        feats += [0.0]

    return np.array(feats, dtype=np.float32)

def get_bond_features(bond):
    bt = bond.GetBondType()
    return np.array([
        int(bt == Chem.rdchem.BondType.SINGLE),
        int(bt == Chem.rdchem.BondType.DOUBLE),
        int(bt == Chem.rdchem.BondType.TRIPLE),
        int(bt == Chem.rdchem.BondType.AROMATIC),
        int(bond.GetIsConjugated()),
        int(bond.IsInRing())
    ], dtype=np.float32)

# ---------- SMILES → graph ----------

def mol_to_graph(smiles: str, label: float = None):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    try:
        rdPartialCharges.ComputeGasteigerCharges(mol)
    except:
        pass

    x = torch.tensor(
        [get_atom_features(a, mol) for a in mol.GetAtoms()],
        dtype=torch.float
    )

    edges, edge_feats = [], []
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bf = get_bond_features(bond)
        edges += [[i,j],[j,i]]
        edge_feats += [bf,bf]

    edge_index = torch.tensor(edges, dtype=torch.long).t().contiguous()
    edge_attr = torch.tensor(edge_feats, dtype=torch.float)

    global_feats = []
    for fn in descriptor_functions:
        try:
            global_feats.append(safe_value(fn(mol)))
        except:
            global_feats.append(0.0)

    data = MoleculeData(
        x=x,
        edge_index=edge_index,
        edge_attr=edge_attr,
        global_features=torch.tensor(global_feats, dtype=torch.float)
    )

    if label is not None:
        data.y = torch.tensor([label], dtype=torch.float)

    return data
