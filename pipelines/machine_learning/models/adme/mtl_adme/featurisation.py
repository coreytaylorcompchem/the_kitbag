import numpy as np
import torch
from torch_geometric.data import Data
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit import DataStructs
from rdkit.Chem import Crippen, rdPartialCharges

from sklearn.preprocessing import StandardScaler, normalize
from sklearn.decomposition import PCA
from joblib import Parallel, delayed

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

from models.adme.utils.compute_fp import compute_morgan_fp

# =========================
# GLOBAL STORAGE (pipeline-safe)
# =========================
GLOBAL_FEATURES = None
FP_FEATURES = None

# =========================
# HELPERS
# =========================

def safe_value(val):
    try:
        if np.isnan(val) or np.isinf(val):
            return 0.0
        return val
    except:
        return 0.0

# =========================
# PRECOMPUTATION
# =========================

def prepare_features(
    df,
    smiles_col="smiles",
    label_col=None,
    global_scaler=None,
    fit_global_scaler=True,
):
    """
    Runs once before featurisation loop.

    If fit_global_scaler=True, fit a new StandardScaler on this dataframe.
    If fit_global_scaler=False, transform descriptors using the supplied pretrained scaler.

    Returns
    -------
    dict
        Contains fitted/used preprocessing objects.
    """
    global GLOBAL_FEATURES, FP_FEATURES

    smiles_list = df[smiles_col].tolist()

    def process_mol(smi):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return np.zeros(28, dtype=np.float32)

        return np.array([
            safe_value(Descriptors.MolWt(mol)),
            safe_value(Descriptors.NumHDonors(mol)),
            safe_value(Descriptors.NumHAcceptors(mol)),
            safe_value(Descriptors.TPSA(mol)),
            safe_value(Descriptors.FractionCSP3(mol)),
            safe_value(Descriptors.HeavyAtomCount(mol)),
            safe_value(Descriptors.RingCount(mol)),
            safe_value(Descriptors.NumRotatableBonds(mol)),
            safe_value(Descriptors.BalabanJ(mol)),
            safe_value(Descriptors.LabuteASA(mol)),
            safe_value(Descriptors.NumAliphaticRings(mol)),
            safe_value(Descriptors.NumAromaticRings(mol)),
            safe_value(Descriptors.NumSaturatedRings(mol)),
            safe_value(Descriptors.NumValenceElectrons(mol)),
            safe_value(Descriptors.MaxPartialCharge(mol)),
            safe_value(Descriptors.MinPartialCharge(mol)),
            safe_value(Descriptors.MaxAbsPartialCharge(mol)),
            safe_value(Descriptors.MinAbsPartialCharge(mol)),
            safe_value(Crippen.MolLogP(mol)),
            safe_value(Crippen.MolMR(mol)),
            safe_value(Descriptors.NumHeteroatoms(mol)),
            safe_value(Descriptors.HeavyAtomMolWt(mol)),
            safe_value(Descriptors.NHOHCount(mol)),
            safe_value(Descriptors.NOCount(mol)),
            safe_value(Descriptors.NumAliphaticHeterocycles(mol)),
            safe_value(Descriptors.NumAromaticHeterocycles(mol)),
            safe_value(Descriptors.NumSaturatedHeterocycles(mol)),
        ], dtype=np.float32)

    global_feats = [process_mol(s) for s in smiles_list]
    global_feats = np.array(global_feats, dtype=np.float32)

    if fit_global_scaler:
        scaler = StandardScaler()
        GLOBAL_FEATURES = scaler.fit_transform(global_feats)
    else:
        if global_scaler is None:
            raise ValueError(
                "fit_global_scaler=False but no global_scaler was supplied."
            )
        scaler = global_scaler
        GLOBAL_FEATURES = scaler.transform(global_feats)

    def process_fp(smi):
        mol = Chem.MolFromSmiles(smi)

        if mol is None:
            return np.zeros(3072, dtype=np.float32)

        def fp_to_numpy(fp, n_bits):
            arr = np.zeros((n_bits,), dtype=np.float32)
            DataStructs.ConvertToNumpyArray(fp, arr)
            return arr

        fp_ecfp = fp_to_numpy(
            AllChem.GetMorganFingerprintAsBitVect(
                mol,
                radius=2,
                nBits=1024,
                useFeatures=False
            ),
            1024
        )

        fp_fcfp = fp_to_numpy(
            AllChem.GetMorganFingerprintAsBitVect(
                mol,
                radius=2,
                nBits=1024,
                useFeatures=True
            ),
            1024
        )

        fp_atompair = fp_to_numpy(
            AllChem.GetHashedAtomPairFingerprintAsBitVect(
                mol,
                nBits=512
            ),
            512
        )

        fp_torsion = fp_to_numpy(
            AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect(
                mol,
                nBits=512
            ),
            512
        )

        return np.concatenate([
            fp_ecfp,
            fp_fcfp,
            fp_atompair,
            fp_torsion
        ])

    fps = [process_fp(s) for s in smiles_list]
    fps = np.array(fps, dtype=np.float32)

    FP_FEATURES = normalize(fps, axis=1)

    return {
        "global_feature_scaler": scaler
    }


# =========================
# GRAPH CONSTRUCTION
# =========================

def mol_to_graph(smiles: str, label=None, idx=None, global_feats=None, fps=None):
    global GLOBAL_FEATURES, FP_FEATURES

    if idx is None:
        raise ValueError("Index required for multitask featurisation")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    Chem.rdPartialCharges.ComputeGasteigerCharges(mol)

    # Node features
    def atom_features(atom):

        try:
            g_charge = float(atom.GetProp('_GasteigerCharge'))
            if np.isnan(g_charge) or np.isinf(g_charge):
                g_charge = 0.0
        except:
            g_charge = 0.0

        return [
            # identity
            atom.GetAtomicNum(),

            # topology
            atom.GetDegree(),
            atom.GetFormalCharge(),
            atom.GetTotalNumHs(),
            atom.GetImplicitValence(),
            atom.GetTotalValence(),
            atom.GetNumRadicalElectrons(),
            atom.GetNumExplicitHs(),

            # electronic
            g_charge,

            # aromaticity / ring
            int(atom.GetIsAromatic()),
            int(atom.IsInRing()),

            int(atom.IsInRingSize(3)),
            int(atom.IsInRingSize(4)),
            int(atom.IsInRingSize(5)),
            int(atom.IsInRingSize(6)),
            int(atom.IsInRingSize(7)),

            # hybridization
            int(atom.GetHybridization() == Chem.rdchem.HybridizationType.SP),
            int(atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2),
            int(atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3),

            # chirality
            int(atom.HasProp('_ChiralityPossible')),

            # medicinal chemistry flags
            int(atom.GetAtomicNum() in [7, 8]),   # hetero
            int(atom.GetAtomicNum() in [9, 17, 35, 53]),  # halogen
            int(Descriptors.NumHDonors(mol) > 0),
            int(Descriptors.NumHAcceptors(mol) > 0),

            int(atom.GetIdx() in donor_matches),
            int(atom.GetIdx() in acceptor_matches),
        ]

    donor_smarts = Chem.MolFromSmarts(
        "[$([N;!H0;v3,v4&+1]),$([O,S;H1;+0]),n&H1&+0]"
    )

    acceptor_smarts = Chem.MolFromSmarts(
        "[$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=[O,N,P,S])]),n&H0&+0]"
    )

    donor_matches = set()
    for match in mol.GetSubstructMatches(donor_smarts):
        donor_matches.update(match)

    acceptor_matches = set()
    for match in mol.GetSubstructMatches(acceptor_smarts):
        acceptor_matches.update(match)

    atom_features_list = [atom_features(atom) for atom in mol.GetAtoms()]
    x = torch.tensor(atom_features_list, dtype=torch.float32)

    # Edges
    edge_index = []
    edge_attr = []

    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()

        edge_index += [[i, j], [j, i]]

        edge_attr += [[
            bond.GetBondTypeAsDouble(),
            int(bond.GetIsConjugated()),
            int(bond.IsInRing()),
            int(bond.GetIsAromatic()),
            int(bond.GetStereo() != Chem.rdchem.BondStereo.STEREONONE),
            int(bond.GetBondType() == Chem.rdchem.BondType.SINGLE),
            int(bond.GetBondType() == Chem.rdchem.BondType.DOUBLE),
            int(bond.GetBondType() == Chem.rdchem.BondType.TRIPLE),
            int(bond.GetBondType() == Chem.rdchem.BondType.AROMATIC),
        ]] * 2

    edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    edge_attr = torch.tensor(edge_attr, dtype=torch.float32)

    data = Data(
        x=x,
        edge_index=edge_index,
        edge_attr=edge_attr,
    )

    # Attach precomputed features
    data.global_features = torch.tensor(GLOBAL_FEATURES[idx]).unsqueeze(0)
    data.fp = torch.tensor(FP_FEATURES[idx]).unsqueeze(0)

    if GLOBAL_FEATURES is None or FP_FEATURES is None:
        raise RuntimeError(
            "GLOBAL_FEATURES / FP_FEATURES not initialised. "
            "Did you forget to call prepare_features()?"
        )

    if label is not None:
        data.y = torch.tensor(label, dtype=torch.float32).view(1, -1)

    return data