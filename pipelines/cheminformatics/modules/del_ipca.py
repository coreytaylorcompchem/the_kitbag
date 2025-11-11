import re, io, contextlib
from collections import Counter

import numpy as np
import pandas as pd
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import ConvertToNumpyArray
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

from sklearn.decomposition import IncrementalPCA

import matplotlib.pyplot as plt
from pathlib import Path

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_task("compute_fingerprints",
               category="DEL Analysis",
               description="Compute Morgan fingerprints for iPCA.")
def compute_fingerprints_task(config, data):
    """
    Compute Morgan fingerprints for molecules and synthons with tqdm progress.
    Handles trailing commas and silences RDKit SMILES parse errors.
    """

    df = data.get("df")
    if df is None or df.empty:
        return {"df": df}

    n_bits = config.get("fingerprint", {}).get("n_bits", 1024)
    radius = config.get("fingerprint", {}).get("radius", 2)

    # Function to clean SMILES robustly (KinDel has issues...)
    def clean_smiles(smiles):
        """
        Extract the first valid-looking SMILES from a string.
        Removes trailing commas, whitespace, quotes.
        """
        if not isinstance(smiles, str):
            return None
        smiles = smiles.strip()
        # Remove quotes
        smiles = smiles.strip('"').strip("'")
        # Remove trailing commas/spaces
        smiles = re.sub(r'[,\s]+$', '', smiles)
        # Extract the first contiguous SMILES substring
        # SMILES allowed characters: alphanumerics, @, =, #, /, \, (), [], +, -, ., :, %
        match = re.match(r'^[A-Za-z0-9@=\#\-/\\\(\)\[\]\+\-\.\:%]+', smiles)
        if match:
            return match.group(0)
        return None

    # Track valid/invalid SMILES
    fp_stats = Counter()

    def mol_fp(smiles: str):
        smiles_clean = clean_smiles(smiles)
        if not smiles_clean:
            fp_stats['invalid'] += 1
            return None
        try:
            # Shut up, RDKit
            with contextlib.redirect_stderr(io.StringIO()):
                mol = Chem.MolFromSmiles(smiles_clean, sanitize=True)
            if mol is None:
                fp_stats['invalid'] += 1
                return None
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits)
            arr = np.zeros((n_bits,), dtype=np.int8)
            ConvertToNumpyArray(fp, arr)
            fp_stats['valid'] += 1
            return arr
        except Exception:
            fp_stats['invalid'] += 1
            return None

    smiles_columns = config.get("compute_fingerprints", {}).get("smiles_columns", {})

    for smiles_col, fp_col in smiles_columns.items():
        if smiles_col not in df.columns:
            logger.warning(f"[compute_fingerprints_task] Missing column '{smiles_col}', skipping.")
            continue
        tqdm.pandas(desc=f"Computing fingerprints for {smiles_col}")
        df[fp_col] = df[smiles_col].progress_apply(mol_fp)
        n_invalid = df[fp_col].isna().sum()
        logger.info(f"[compute_fingerprints_task] {smiles_col}: {n_invalid:,} invalid SMILES skipped.")

    logger.info(f"[compute_fingerprints_task] Fingerprint summary: {fp_stats}")

    # Only drop rows for fingerprint columns that exist
    subset = [fp_col for fp_col in ["fp_a", "fp_b", "fp_c", "fp_final"] if fp_col in df.columns]
    if subset:
        df = df.dropna(subset=subset).reset_index(drop=True)

    # Concatenate fp cols columns
    df["fingerprints"] = df.apply(
        lambda row: np.concatenate([row[c] for c in subset if row[c] is not None]),
        axis=1
    )

    df = df.head(1000)

    logger.info(f"[compute_fingerprints_task] Computed fingerprints for {len(df)} molecules.")
    return {"df": df}


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

    # === Prototype slicing ===
    proto_config = config.get("compute_fingerprints", {}).get("prototyping", {})
    if proto_config:
        n_rows = proto_config.get("n_rows", None)
        random_sample = proto_config.get("random_sample", False)
        if n_rows and n_rows < len(df):
            if random_sample:
                df = df.sample(n=n_rows, random_state=42).reset_index(drop=True)
            else:
                df = df.head(n_rows).copy()
            logger.info(f"[compute_fingerprints_task] Using prototype subset: {len(df)} rows")

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