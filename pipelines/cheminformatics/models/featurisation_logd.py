import itertools
import os
import pathlib
from pathlib import Path

from multiprocessing import Pool

import pandas as pd
import numpy as np

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, Descriptors3D, rdMolDescriptors, rdPartialCharges
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

import torch
import torch.nn.functional as F
import torch.nn as nn

# from torch_geometric
from torch_geometric.data import Data

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

def get_morgan_fingerprint(mol, radius=2, nBits=1024):
    generator = GetMorganGenerator(radius=radius, fpSize=nBits)
    fp = generator.GetFingerprint(mol)
    arr = np.zeros((nBits,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr.astype(np.float32)

def get_3d_descriptors(mol):
    mol = Chem.AddHs(mol)

    try:
        if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
            return np.zeros(6, dtype=np.float32)

        # Check for UFF optimization support
        try:
            AllChem.UFFOptimizeMolecule(mol)
        except:
            return np.zeros(6, dtype=np.float32)

        descs = [
            Descriptors3D.Asphericity(mol),
            Descriptors3D.Eccentricity(mol),
            Descriptors3D.InertialShapeFactor(mol),
            Descriptors3D.SpherocityIndex(mol),
            rdMolDescriptors.CalcVolume(mol),
            rdMolDescriptors.CalcLabuteASA(mol)
        ]

        descs = [0.0 if np.isnan(d) or np.isinf(d) else d for d in descs]
        return np.array(descs, dtype=np.float32)

    except Exception as e:
        return np.zeros(6, dtype=np.float32)

# We can get reasonable results (R^2 ~ 0.5) with simple atomic features.
# But as we are trying to predict a more complex phenomena, we can benefit from RdKit 
# descriptors to encode hybridisation, partial charges, etc.
# LogD prediction should improve from here.

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

# Precompute a dict of electronegativities (Pauling scale)
electronegativity_dict = {
    1: 2.20, 6: 2.55, 7: 3.04, 8: 3.44, 9: 3.98, 15: 2.19,
    16: 2.58, 17: 3.16, 35: 2.96, 53: 2.66,
}

# List of metals by atomic number (common metals)
metals = {
    3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 30, 37, 38, 39, 40, 41, 42,
    43, 44, 45, 46, 47, 48, 49, 55, 56, 57, 72,
    73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83,
    87, 88, 89, 104, 105, 106, 107, 108, 109, 110,
    111, 112, 113, 114, 115, 116
}

def one_hot_encoding(x, allowable_set):
    if x not in allowable_set:
        x = allowable_set[-1]  # fallback to "unknown"
    return [int(x == s) for s in allowable_set]

def get_atom_features(atom, mol=None):
    pt = Chem.GetPeriodicTable()
    features = []

    # Atom symbol (one-hot + unknown)
    features += one_hot_encoding(atom.GetSymbol(),
                                ['C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I', 'H',
                                 'B', 'Si', 'Se', 'As', 'Al', 'Zn', 'Cu', 'Ni', 'Fe', 'other'])

    # Hybridization (one-hot)
    features += one_hot_encoding(atom.GetHybridization(),
                                [Chem.rdchem.HybridizationType.SP,
                                 Chem.rdchem.HybridizationType.SP2,
                                 Chem.rdchem.HybridizationType.SP3,
                                 Chem.rdchem.HybridizationType.SP3D,
                                 Chem.rdchem.HybridizationType.SP3D2,
                                 Chem.rdchem.HybridizationType.UNSPECIFIED])

    # Integer-valued features (can be scaled later if needed)
    features += [
        atom.GetDegree(),               # n_bonded neighbors
        atom.GetFormalCharge(),
        atom.GetNumRadicalElectrons(),
        int(atom.GetIsAromatic())
    ]

    atomic_num = atom.GetAtomicNum()

    # Physical properties from periodic table (normalised to 0-1)
    # We could prpbably play with this, if we want to.
    features += [
        pt.GetAtomicWeight(atomic_num) / 200,                 # ~1–200 → 0–1
        pt.GetRvdw(atomic_num) / 2.5,                         # ~0.5–2.5 → 0–1
        electronegativity_dict.get(atomic_num, 0.0) / 4.0,    # ~0–4 → 0–1
        pt.GetNOuterElecs(atomic_num) / 8.0,                  # Valence electrons, ~0–8
        1.0 if atomic_num in metals else 0.0                  # binary
    ]


    # Partial charge (Gasteiger)
    # Something else we could play with (e.g. AM1-BCC) but this adds cost.
    if mol is not None:
        try:
            charge = float(atom.GetProp('_GasteigerCharge'))
            if np.isnan(charge) or np.isinf(charge):
                charge = 0.0
        except:
            charge = 0.0
        features += [charge]
    else:
        features += [0.0]

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
    ]
    return np.array(features, dtype=np.float32)

class MoleculeData(Data):
    def __init__(self, x, edge_index, edge_attr, global_features):
        super().__init__(x=x, edge_index=edge_index, edge_attr=edge_attr)
        self.global_features = global_features

    def __cat_dim__(self, key, value, *args, **kwargs):
        if key == 'global_features':
            return None
        return super().__cat_dim__(key, value)

    def to(self, device):
        out = super().to(device)  # Let torch_geometric handle all attributes
        out.global_features = self.global_features.to(device)
        return out

def mol_to_graph(smiles: str, label: float = None, radius=2, nBits=1024):
    print(f"[mol_to_graph] Processing SMILES: {smiles}")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"[mol_to_graph] RDKit failed to parse: {smiles}")
        return None

    try:
        rdPartialCharges.ComputeGasteigerCharges(mol)
        print("[mol_to_graph] Gasteiger charges computed.")
    except Exception as e:
        print(f"[mol_to_graph] Warning: Gasteiger charges failed for {smiles}: {e}")

    try:
        atom_features = [get_atom_features(atom, mol) for atom in mol.GetAtoms()]
        x = torch.tensor(atom_features, dtype=torch.float)
        print(f"[mol_to_graph] Atom features shape: {x.shape}")
    except Exception as e:
        print(f"[mol_to_graph] Failed to get atom features: {e}")
        return None

    try:
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
        print(f"[mol_to_graph] Edge index shape: {edge_index.shape}, edge_attr shape: {edge_attr.shape}")
    except Exception as e:
        print(f"[mol_to_graph] Failed to get edge features: {e}")
        return None

    try:
        global_features = []

        for func in descriptor_functions:
            try:
                val = func(mol)
                if np.isnan(val) or np.isinf(val):
                    val = 0.0
            except:
                val = 0.0
            global_features.append(val)

        fp = get_morgan_fingerprint(mol, radius=radius, nBits=nBits)
        global_features += fp.tolist()

        desc3d = get_3d_descriptors(mol)
        global_features += desc3d.tolist()

        global_features = torch.tensor(global_features, dtype=torch.float32)
        print(f"[mol_to_graph] Global features shape: {global_features.shape}")

        if not isinstance(global_features, torch.Tensor):
            print(f"[mol_to_graph] WARNING: global_features is not a tensor, type={type(global_features)}")
        data = MoleculeData(x=x, edge_index=edge_index, edge_attr=edge_attr, global_features=global_features)

        if label is not None:
            data.y = torch.tensor([label], dtype=torch.float)

        print(f"[mol_to_graph] Successfully created graph for: {smiles}")
        return data

    except Exception as e:
        print(f"[mol_to_graph] Failed to create MoleculeData object for {smiles}: {e}")
        return None
