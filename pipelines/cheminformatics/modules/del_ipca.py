import numpy as np
import pandas as pd
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import ConvertToNumpyArray

from sklearn.decomposition import IncrementalPCA

import matplotlib.pyplot as plt
from pathlib import Path

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_task("compute_fingerprints",
               category="DEL Analysis",
               description="Compute fps for iPCA.")
def compute_fingerprints_task(config, data):
    """
    Compute Morgan fingerprints for molecules and synthons with tqdm progress.
    """
    df = data.get("df")
    if df is None or df.empty:
        return {"df": df}

    n_bits = config.get("fingerprint", {}).get("n_bits", 1024)
    radius = config.get("fingerprint", {}).get("radius", 2)

    def mol_fp(smiles: str):
        try:
            smiles_clean = smiles.strip().rstrip(',')
            mol = Chem.MolFromSmiles(smiles_clean, sanitize=True)
            if mol is None:
                return None
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)
            arr = np.zeros((1024,), dtype=np.int8)
            ConvertToNumpyArray(fp, arr)
            return arr
        except Exception:
            # Any parsing errors are silently returned as None
            return None

    tqdm.pandas(desc="Computing fingerprints")  # wrap apply with progress bar

    for col in ["smiles_a", "smiles_b", "smiles_c", "smiles"]:
        df[f"fp_{col.split('_')[-1]}"] = df[col].progress_apply(mol_fp)

    # Drop failed rows
    df = df.dropna(subset=["fp_a", "fp_b", "fp_c", "fp_final"]).reset_index(drop=True)

    # Concatenate into single feature vector
    df["fingerprints"] = df.apply(lambda row: np.concatenate([row["fp_a"], row["fp_b"], row["fp_c"], row["fp_final"]]), axis=1)

    logger.info(f"[compute_fingerprints_task] Computed fingerprints for {len(df)} molecules.")
    return {"df": df}

# -------------------------------
@register_task("incremental_pca",
               category="DEL Analysis",
               description="Initial dimensionality reduction to produce reduced dataset.")
def incremental_pca_task(config, data):
    """
    Run incremental PCA on fingerprint columns with batch progress.
    """
    df = data.get("df")
    if df is None or df.empty:
        return {"df": df}

    X = np.stack(df["fingerprints"].values)

    ipca = config.get("_ipca_instance")
    if ipca is None:
        n_components = config.get("incremental_pca", {}).get("n_components", 50)
        batch_size = config.get("incremental_pca", {}).get("batch_size", 50000)
        ipca = IncrementalPCA(n_components=n_components, batch_size=batch_size)
        config["_ipca_instance"] = ipca

    # Partial fit
    logger.info(f"[incremental_pca_task] Performing partial_fit on {X.shape[0]} samples.")
    ipca.partial_fit(X)

    # Transform chunk
    X_reduced = ipca.transform(X)
    for i in range(X_reduced.shape[1]):
        df[f"PC{i+1}"] = X_reduced[:, i]

    return {"df": df}