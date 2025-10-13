import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem import Descriptors, rdPartialCharges
from torch_geometric.data import Data

# ---- Custom Data class ----

class MoleculeData(Data):
    def __inc__(self, key, value, *args, **kwargs):
        return super().__inc__(key, value)

    def __cat_dim__(self, key, value, *args, **kwargs):
        if key == 'global_features':
            return None  # do not concatenate along any dimension
        return super().__cat_dim__(key, value)

# ---- Descriptor and feature helpers ----

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
    1: 2.20, 6: 2.55, 7: 3.04, 8: 3.44, 9: 3.98, 15: 2.19,
    16: 2.58, 17: 3.16, 35: 2.96, 53: 2.66,
}

polarizability_dict = {
    1: 0.666, 6: 1.76, 7: 1.1, 8: 0.802, 9: 0.557, 15: 3.63,
    16: 3.6, 17: 2.18,
}

metals = {
    3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 30, 37, 38, 39, 40, 41, 42,
    43, 44, 45, 46, 47, 48, 49, 55, 56, 57, 72,
    73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83,
    87, 88, 89, 104, 105, 106, 107, 108, 109, 110,
    111, 112, 113, 114, 115, 116
}

def safe_value(x):
    if x is None or np.isnan(x) or np.isinf(x):
        return 0.0
    return x

def one_hot_encoding(x, allowable_set):
    if x not in allowable_set:
        x = allowable_set[-1]
    return [int(x == s) for s in allowable_set]

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
        atom.GetFormalCharge() / 3.0,
        atom.GetNumRadicalElectrons(),
        int(atom.GetIsAromatic()),
        int(atom.IsInRing()),
    ]

    ring_sizes = [3, 4, 5, 6, 7, 8]
    for size in ring_sizes:
        features.append(int(atom.IsInRingSize(size)))

    atomic_num = atom.GetAtomicNum()
    features += [
        pt.GetAtomicWeight(atomic_num) / 200,
        pt.GetRvdw(atomic_num) / 2.5,
        electronegativity_dict.get(atomic_num, 0.0) / 4.0,
        polarizability_dict.get(atomic_num, 1.0) / 5.0,
        pt.GetNOuterElecs(atomic_num) / 8.0,
        1.0 if atomic_num in metals else 0.0
    ]

    if mol is not None:
        try:
            charge = float(atom.GetProp('_GasteigerCharge'))
            charge = safe_value(charge)
        except:
            charge = 0.0
        features.append(charge)
    else:
        features.append(0.0)

    return np.array(features, dtype=np.float32)

def get_bond_features(bond):
    bond_type = bond.GetBondType()
    features = [
        int(bond_type == Chem.rdchem.BondType.SINGLE),
        int(bond_type == Chem.rdchem.BondType.DOUBLE),
        int(bond_type == Chem.rdchem.BondType.TRIPLE),
        int(bond_type == Chem.rdchem.BondType.AROMATIC),
        int(bond.GetIsConjugated()),
        int(bond.IsInRing()),
        int(bond.GetStereo() == Chem.rdchem.BondStereo.STEREOZ),
        int(bond.GetStereo() == Chem.rdchem.BondStereo.STEREOE),
        int(bond.GetStereo() == Chem.rdchem.BondStereo.STEREOANY),
    ]
    return np.array(features, dtype=np.float32)

metabolic_smarts = {
    'aromatic_amine': '[a][NH2,NH,N]',
    'phenol': 'c[OH]',
    'ester': 'C(=O)O',
    'ether': 'C-O-C',
    'thioether': 'C-S-C',
    'nitro': '[N+](=O)[O-]',
    'heteroaromatic': '[a;!#6]',
}

compiled_smarts = {k: Chem.MolFromSmarts(v) for k, v in metabolic_smarts.items()}

def get_metabolic_flags(mol):
    return [int(mol.HasSubstructMatch(pattern)) for pattern in compiled_smarts.values()]

def mol_to_graph(smiles: str, label=None):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    try:
        rdPartialCharges.ComputeGasteigerCharges(mol)
    except Exception:
        pass

    atom_features = [get_atom_features(atom, mol) for atom in mol.GetAtoms()]
    x = torch.tensor(atom_features, dtype=torch.float)

    edge_index = []
    edge_attr = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        bond_feat = get_bond_features(bond)
        edge_index += [[i, j], [j, i]]
        edge_attr += [bond_feat, bond_feat]

    edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    edge_attr = torch.tensor(edge_attr, dtype=torch.float)

    global_features = []
    for func in descriptor_functions:
        try:
            val = func(mol)
            val = safe_value(val)
        except Exception:
            val = 0.0
        global_features.append(val)

    global_features += get_metabolic_flags(mol)
    global_features = torch.tensor(global_features, dtype=torch.float32)

    data = MoleculeData(x=x, edge_index=edge_index, edge_attr=edge_attr, global_features=global_features)

    if label is not None:
        data.y = torch.tensor([label], dtype=torch.float)

    return data