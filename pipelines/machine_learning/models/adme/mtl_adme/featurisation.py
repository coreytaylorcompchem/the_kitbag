import numpy as np
import torch
from torch_geometric.data import Data
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit import DataStructs

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

def prepare_features(df, smiles_col="smiles", label_col=None):
    """
    Runs ONCE before featurisation loop.
    Stores results in module-level globals.
    """
    global GLOBAL_FEATURES, FP_FEATURES

    smiles_list = df[smiles_col].tolist()

    # ---------------------------
    # GLOBAL DESCRIPTORS
    # ---------------------------
    def process_mol(smi):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return np.zeros(16, dtype=np.float32)

        return np.array([
            safe_value(Descriptors.MolWt(mol)),
            safe_value(Descriptors.MolLogP(mol)),
            safe_value(Descriptors.NumHDonors(mol)),
            safe_value(Descriptors.NumHAcceptors(mol)),
            safe_value(Descriptors.TPSA(mol)),
        ], dtype=np.float32)

    global_feats = [process_mol(s) for s in smiles_list]
    global_feats = np.array(global_feats, dtype=np.float32)

    scaler = StandardScaler()
    GLOBAL_FEATURES = scaler.fit_transform(global_feats)

    # ---------------------------
    # FINGERPRINTS + PCA
    # ---------------------------
    def process_fp(smi):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return np.zeros(1024, dtype=np.float32)

        fp_ecfp = compute_morgan_fp(mol, 512)
        fp_torsion = np.array(
            AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect(mol, nBits=512)
        )
        return np.concatenate([fp_ecfp, fp_torsion])

    fps = [process_fp(s) for s in smiles_list]
    fps = np.array(fps, dtype=np.float32)

    # Remove PCA of fps for now
    
    # pca = PCA(n_components=512)
    # FP_FEATURES = normalize(pca.fit_transform(fps), axis=1)

    FP_FEATURES = normalize(fps, axis=1)


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

    # Node features
    def atom_features(atom):
        return [
            atom.GetAtomicNum(),
            atom.GetDegree(),
            atom.GetFormalCharge(),
            int(atom.GetIsAromatic()),
            atom.GetHybridization().real,
            atom.GetTotalNumHs(),
            atom.GetImplicitValence(),
            int(atom.IsInRing()),
        ]

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
            int(bond.GetBondTypeAsDouble()),
            int(bond.GetIsConjugated()),
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