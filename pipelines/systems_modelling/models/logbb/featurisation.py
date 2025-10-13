import numpy as np
import torch
from torch_geometric.data import Data, Batch
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, rdPartialCharges
from rdkit.Chem.rdchem import ChiralType

from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)


# ==========================================================
# Custom Data structure
# ==========================================================
class MoleculeData(Data):
    def __inc__(self, key, value, *args, **kwargs):
        return super().__inc__(key, value)

    def __cat_dim__(self, key, value, *args, **kwargs):
        if key == "global_features":
            return None
        return super().__cat_dim__(key, value)


# ==========================================================
# Helpers and constants
# ==========================================================
def safe_value(x):
    try:
        if np.isnan(x) or np.isinf(x):
            return 0.0
        return float(x)
    except:
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
    3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,
    29, 30, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
    49, 55, 56, 57, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81,
    82, 83, 87, 88, 89, 104, 105, 106, 107, 108, 109, 110,
    111, 112, 113, 114, 115, 116
}


# ==========================================================
# Atom features
# ==========================================================

def get_atom_features(atom, mol=None):
    pt = Chem.GetPeriodicTable()
    features = []

    # --- Basic atom type & hybridization
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

    # --- Local atomic environment
    features += [
        atom.GetDegree(),
        atom.GetFormalCharge(),
        atom.GetNumRadicalElectrons(),
        int(atom.GetIsAromatic())
    ]

    # --- Aromatic neighbors
    neigh_aromatic = sum(1 for n in atom.GetNeighbors() if n.GetIsAromatic())
    features += [neigh_aromatic / max(1, atom.GetDegree())]

    # --- Ring membership (sizes 3–6)
    ring_info = mol.GetRingInfo() if mol is not None else None
    if ring_info:
        for size in [3, 4, 5, 6]:
            features += [int(ring_info.IsAtomInRingOfSize(atom.GetIdx(), size))]
    else:
        features += [0, 0, 0, 0]

    # --- Physical/chemical constants (normalized)
    atomic_num = atom.GetAtomicNum()
    features += [
        pt.GetAtomicWeight(atomic_num) / 150,                 # normalized mass
        pt.GetRvdw(atomic_num) / 2.5,
        electronegativity_dict.get(atomic_num, 2.5) / 4.0,
        pt.GetNOuterElecs(atomic_num) / 8.0,
        1.0 if atomic_num in metals else 0.0,
        atom.GetTotalNumHs(includeNeighbors=True) / 4.0
    ]

    # --- Implicit hydrogens count (new)
    features.append(atom.GetNumImplicitHs())

    # --- Chirality one-hot encoding (new)
    chiral_tag = atom.GetChiralTag()
    chiral_one_hot = one_hot_encoding(chiral_tag, [
        ChiralType.CHI_UNSPECIFIED,
        ChiralType.CHI_TETRAHEDRAL_CW,
        ChiralType.CHI_TETRAHEDRAL_CCW,
        ChiralType.CHI_OTHER
    ])
    features += chiral_one_hot

    # --- Gasteiger charge (approximate partial charge)
    if mol is not None:
        try:
            charge = float(atom.GetProp('_GasteigerCharge'))
            charge = safe_value(charge)
        except:
            charge = 0.0
        features += [charge]
    else:
        features += [0.0]

    return np.array(features, dtype=np.float32)


# ==========================================================
# Bond features
# ==========================================================

def get_bond_features(bond, mol=None):
    i, j = bond.GetBeginAtom(), bond.GetEndAtom()
    bond_type = bond.GetBondType()

    # --- Basic bond types & conjugation
    features = [
        int(bond_type == Chem.rdchem.BondType.SINGLE),
        int(bond_type == Chem.rdchem.BondType.DOUBLE),
        int(bond_type == Chem.rdchem.BondType.TRIPLE),
        int(bond_type == Chem.rdchem.BondType.AROMATIC),
        int(bond.GetIsConjugated()),
        int(bond.IsInRing()),
    ]

    # --- Stereo and bond order
    features += [
        int(bond.GetStereo() != Chem.rdchem.BondStereo.STEREONONE),
        bond.GetBondTypeAsDouble() / 3.0
    ]

    # --- Electronegativity difference
    en_diff = abs(electronegativity_dict.get(i.GetAtomicNum(), 2.5) -
                  electronegativity_dict.get(j.GetAtomicNum(), 2.5))
    features += [en_diff / 4.0]

    # --- Is bond both aromatic and conjugated (new)
    features.append(int(bond.GetIsConjugated() and bond_type == Chem.rdchem.BondType.AROMATIC))

    # --- Sum of degrees of connected atoms (new)
    deg_sum = i.GetDegree() + j.GetDegree()
    features.append(deg_sum / 20)  # normalized by 20 as an example max degree

    # --- Bond ring size (new)
    if mol is not None:
        ring_info = mol.GetRingInfo()
        bond_rings = ring_info.BondRings()
        rings_containing_bond = [ring for ring in bond_rings if bond.GetIdx() in ring]
        if rings_containing_bond:
            smallest_ring_size = min(len(ring) for ring in rings_containing_bond)
        else:
            smallest_ring_size = 0
        features.append(smallest_ring_size / 6)  # normalize by max ring size 6
    else:
        features.append(0)

    # --- Bond length (if 3D coords available, new)
    if mol is not None and mol.GetNumConformers() > 0:
        pos = mol.GetConformer().GetPositions()
        bond_length = np.linalg.norm(pos[bond.GetBeginAtomIdx()] - pos[bond.GetEndAtomIdx()])
        features.append(bond_length / 1.5)  # normalize by ~1.5 Å typical bond length
    else:
        features.append(0)

    return np.array(features, dtype=np.float32)



# ==========================================================
# Global molecular descriptors
# ==========================================================

descriptor_functions = [
    Descriptors.MolWt,
    Descriptors.MolLogP,
    Descriptors.NumHDonors,
    Descriptors.NumHAcceptors,
    Descriptors.TPSA,
    Descriptors.FractionCSP3,
    Descriptors.HeavyAtomCount,
    Descriptors.NumRotatableBonds,
    Descriptors.RingCount,
    rdMolDescriptors.CalcNumAmideBonds,
    rdMolDescriptors.CalcLabuteASA,
    rdMolDescriptors.CalcKappa1,
    rdMolDescriptors.CalcChi0n,
    rdMolDescriptors.CalcChi1n,
    lambda mol: rdMolDescriptors.CalcCrippenDescriptors(mol)[0],  # logP
    lambda mol: rdMolDescriptors.CalcCrippenDescriptors(mol)[1],  # MR
    rdMolDescriptors.CalcLabuteASA,
    # rdMolDescriptors.CalcChi0n, rdMolDescriptors.CalcChi1n,
    # rdMolDescriptors.CalcChi0v, rdMolDescriptors.CalcChi1v
]


# ==========================================================
# Molecule → Graph
# ==========================================================
def mol_to_graph(smiles: str, label: float = None):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # Ensure the molecule is sanitised so ring info and valences are computed
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        logger.warning(f"Molecule sanitisation failed for {smiles}: {e}")

    try:
        rdPartialCharges.ComputeGasteigerCharges(mol)
    except Exception:
        pass

    x = torch.tensor([get_atom_features(a, mol) for a in mol.GetAtoms()], dtype=torch.float)
    edge_index, edge_attr = [], []
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bfeat = get_bond_features(bond, mol)
        edge_index += [[i, j], [j, i]]
        edge_attr += [bfeat, bfeat]

    edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    edge_attr = torch.tensor(edge_attr, dtype=torch.float)

    global_features = []
    for func in descriptor_functions:
        try:
            val = safe_value(func(mol))
        except Exception:
            val = 0.0
        global_features.append(val)
    global_features = torch.tensor(global_features, dtype=torch.float32)

    data = MoleculeData(
        x=x, edge_index=edge_index, edge_attr=edge_attr, global_features=global_features
    )

    if label is not None:
        data.y = torch.tensor([label], dtype=torch.float)

    data.smiles = smiles
    return data


# ==========================================================
# Loader utilities
# ==========================================================
def compute_global_features(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return torch.zeros(len(descriptor_functions), dtype=torch.float)
    feats = [safe_value(fn(mol)) for fn in descriptor_functions]
    return torch.tensor(feats, dtype=torch.float)


def custom_collate(batch_list):
    return Batch.from_data_list(batch_list)
