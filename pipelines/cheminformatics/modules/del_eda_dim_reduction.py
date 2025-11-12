import logging
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import re
import contextlib
import io

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

# -----------------------------
# Helpers
# -----------------------------

def clean_smiles(smiles: str) -> str | None:
    """
    Robust SMILES cleaning:
    - Strip whitespace and quotes
    - Remove trailing commas/spaces
    - Extract first contiguous valid SMILES substring
    """
    if not isinstance(smiles, str):
        return None
    smiles = smiles.strip().strip('"').strip("'")
    smiles = re.sub(r'[,\s]+$', '', smiles)
    match = re.match(r'^[A-Za-z0-9@=\#\-/\\\(\)\[\]\+\-\.\:%]+', smiles)
    if match:
        return match.group(0)
    return None

def compute_fingerprint(smiles, radius=2, n_bits=1024):
    smiles_clean = clean_smiles(smiles)
    if not smiles_clean:
        return None
    try:
        with contextlib.redirect_stderr(io.StringIO()):  # silence RDKit warnings
            mol = Chem.MolFromSmiles(smiles_clean, sanitize=True)
        if mol is None:
            return None
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
        arr = np.zeros((n_bits,), dtype=np.int8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr
    except Exception:
        return None

def get_fingerprints(df):
    df['fp_a'] = df['smiles_a'].apply(compute_fingerprint)
    df['fp_b'] = df['smiles_b'].apply(compute_fingerprint)
    df['fp_c'] = df['smiles_c'].apply(compute_fingerprint)
    df['fp_final'] = df['smiles'].apply(compute_fingerprint)
    df = df.dropna(subset=['fp_a', 'fp_b', 'fp_c', 'fp_final']).reset_index(drop=True)
    return df

def plot_2d_projection(X, color=None, title="Dimensionality Reduction",
                       x_label="Component 1", y_label="Component 2",
                       cmap='coolwarm', output_file=None):
    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(X[:, 0], X[:, 1], c=color, cmap=cmap,
                          alpha=0.8, edgecolors='w', s=80)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if color is not None:
        plt.colorbar(scatter, label='Color scale')
    if output_file:
        plt.savefig(output_file, bbox_inches='tight', dpi=300)
        plt.close()
        logger.info(f"Saved plot: {output_file}")
    else:
        plt.show()

# -----------------------------
# Main Task: Dimensionality Reduction
# -----------------------------

@register_task("dimensionality_reduction_analyses",
               category="DEL Analysis",
               description="Run PCA, t-SNE, and UMAP on DEL fingerprints and plot enrichment maps.")
def dimensionality_reduction_analyses(config, data=None):
    input_file = config.get("input_file")
    output_dir = Path(config.get("dimensionality_reduction_analyses", {}).get("output_dir", "outputs/dimensionality_reduction"))
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"[DEL] Reading Parquet file: {input_file}")
    df = pd.read_parquet(input_file)

    # Compute fingerprints robustly
    df = get_fingerprints(df)

    # Concatenate bit vectors for combined feature representation
    fingerprints = np.array([
        np.concatenate([row['fp_a'], row['fp_b'], row['fp_c'], row['fp_final']])
        for _, row in df.iterrows()
    ])

    # Normalize features
    scaler = StandardScaler()
    fingerprints_scaled = scaler.fit_transform(fingerprints)

    # Compute enrichment
    df['target_mean'] = df[['seq_target_1', 'seq_target_2', 'seq_target_3']].mean(axis=1)
    df['matrix_mean'] = df[['seq_matrix_1', 'seq_matrix_2', 'seq_matrix_3']].mean(axis=1)
    df['log_enrichment'] = np.log10((df['target_mean'] + 1) / (df['matrix_mean'] + 1))

    # --- PCA ---
    pca_params = config["dimensionality_reduction_analyses"].get("pca", {})
    n_pca = pca_params.get("n_components", 2)
    pca = PCA(n_components=n_pca)
    pca_result = pca.fit_transform(fingerprints_scaled)
    plot_2d_projection(
        pca_result, color=df['log_enrichment'],
        title="PCA colored by log10(target/control)",
        output_file=output_dir / "pca_projection.png"
    )

    # --- t-SNE ---
    tsne_params = config["dimensionality_reduction_analyses"].get("tsne", {})
    tsne = TSNE(
        n_components=tsne_params.get("n_components", 2),
        perplexity=tsne_params.get("perplexity", 30),
        max_iter=tsne_params.get("max_iter", 300)
    )
    tsne_result = tsne.fit_transform(fingerprints_scaled)
    plot_2d_projection(
        tsne_result, color=df['log_enrichment'],
        title="t-SNE colored by log10(target/control)",
        output_file=output_dir / "tsne_projection.png"
    )

    # --- UMAP ---
    umap_params = config["dimensionality_reduction_analyses"].get("umap", {})
    umap_model = umap.UMAP(
        n_components=umap_params.get("n_components", 2),
        metric=umap_params.get("metric", "jaccard"),
        n_neighbors=umap_params.get("n_neighbors", 30),
        min_dist=umap_params.get("min_dist", 0.1),
        random_state=umap_params.get("random_state", 42)
    )
    umap_result = umap_model.fit_transform(fingerprints_scaled)
    plot_2d_projection(
        umap_result, color=df['log_enrichment'],
        title="UMAP (Jaccard) colored by log10(target/control)",
        output_file=output_dir / "umap_projection.png"
    )

    # Save reduced coordinates
    reduced_df = df.copy()
    reduced_df[['pca_1', 'pca_2']] = pca_result
    reduced_df[['tsne_1', 'tsne_2']] = tsne_result
    reduced_df[['umap_1', 'umap_2']] = umap_result

    output_file = output_dir / "dimensionality_reduction_results.parquet"
    reduced_df.to_parquet(output_file)
    logger.info(f"Saved dimensionality reduction data: {output_file}")

    return {
        "pca": str(output_dir / "pca_projection.png"),
        "tsne": str(output_dir / "tsne_projection.png"),
        "umap": str(output_dir / "umap_projection.png"),
        "reduced_data": str(output_file)
    }
