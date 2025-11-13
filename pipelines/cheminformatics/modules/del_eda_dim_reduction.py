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
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
import umap

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

# Helpers

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

def plot_2d_projection(X, color=None, color_label=None, title="Dimensionality Reduction",
                       x_label="Component 1", y_label="Component 2",
                       cmap='coolwarm', output_file=None):
    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(
        X[:, 0], X[:, 1],
        c=color, cmap=cmap, alpha=0.8,
        edgecolors='w', s=60, linewidths=0.5
    )
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    if color is not None:
        label = color_label if color_label else "Value"
        cbar = plt.colorbar(scatter)
        cbar.set_label(label, rotation=270, labelpad=15)

    if output_file:
        plt.savefig(output_file, bbox_inches='tight', dpi=300)
        plt.close()
        logger.info(f"Saved plot: {output_file}")
    else:
        plt.show()

@register_task("dimensionality_reduction_analyses",
               category="DEL Analysis",
               description="Run PCA, t-SNE, and UMAP on DEL fingerprints, "
                           "analyze enrichment patterns, and produce summary plots.")
def dimensionality_reduction_analyses(config, data=None):
    input_file = config.get("input_file")
    output_dir = Path(config.get("dimensionality_reduction_analyses", {}).get("output_dir", "outputs/dimensionality_reduction"))
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"[DEL] Reading Parquet file: {input_file}")
    df = pd.read_parquet(input_file)

    df = get_fingerprints(df)

    # --- Fingerprint weighting --------------------------------------------------
    weights = config["dimensionality_reduction_analyses"].get("fp_weights", [1, 1, 1, 1])
    if len(weights) != 4:
        logger.warning("fp_weights must have 4 elements (a,b,c,final). Using equal weights.")
        weights = [1, 1, 1, 1]

    fingerprints = np.array([
        np.concatenate([
            weights[0] * row['fp_a'],
            weights[1] * row['fp_b'],
            weights[2] * row['fp_c'],
            weights[3] * row['fp_final']
        ])
        for _, row in df.iterrows()
    ])

    # --- Normalise features -----------------------------------------------------
    scaler = StandardScaler()
    fingerprints_scaled = scaler.fit_transform(fingerprints)

    # --- Compute enrichment -----------------------------------------------------
    df['target_mean'] = df[['seq_target_1', 'seq_target_2', 'seq_target_3']].mean(axis=1)
    df['matrix_mean'] = df[['seq_matrix_1', 'seq_matrix_2', 'seq_matrix_3']].mean(axis=1)
    df['log_enrichment'] = np.log10((df['target_mean'] + 1) / (df['matrix_mean'] + 1))

    # --- PCA --------------------------------------------------------------------
    pca_params = config["dimensionality_reduction_analyses"].get("pca", {})
    n_pca = pca_params.get("n_components", 2)
    pca = PCA(n_components=n_pca)
    pca_result = pca.fit_transform(fingerprints_scaled)

    plot_2d_projection(
        pca_result, color=df['log_enrichment'],
        color_label="log10(target/control)",
        title="PCA colored by enrichment factor",
        output_file=output_dir / "pca_projection.png"
    )

    # --- PCA–Enrichment correlation ---------------------------------------------
    if config["dimensionality_reduction_analyses"].get("pca_correlation", False):
        corrs = [
            np.corrcoef(pca_result[:, i], df["log_enrichment"])[0, 1]
            for i in range(pca_result.shape[1])
        ]
        corr_df = pd.DataFrame({"PC": range(1, len(corrs) + 1),
                                "corr_with_enrichment": corrs})
        corr_df.to_csv(output_dir / "pca_enrichment_correlation.csv", index=False)
        logger.info("Saved PCA–enrichment correlation file.")

    # --- t-SNE ------------------------------------------------------------------
    tsne_params = config["dimensionality_reduction_analyses"].get("tsne", {})
    tsne = TSNE(
        n_components=tsne_params.get("n_components", 2),
        perplexity=tsne_params.get("perplexity", 30),
        max_iter=tsne_params.get("max_iter", 500),
        init=tsne_params.get("init", "pca"),
        learning_rate=tsne_params.get("learning_rate", "auto")
    )
    tsne_result = tsne.fit_transform(fingerprints_scaled)
    plot_2d_projection(
        tsne_result, color=df['log_enrichment'],
        color_label="log10(target/control)",
        title="t-SNE colored by enrichment factor",
        output_file=output_dir / "tsne_projection.png"
    )

    # --- UMAP -------------------------------------------------------------------
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
        color_label="log10(target/control)",
        title="UMAP (Jaccard) colored by enrichment factor",
        output_file=output_dir / "umap_projection.png"
    )

    # --- Property-based coloring ------------------------------------------------
    plot_props = config["dimensionality_reduction_analyses"].get("plot_properties", [])
    if plot_props:
        for prop in plot_props:
            if prop == "MolWt":
                df[prop] = df["smiles"].apply(lambda s: Descriptors.MolWt(Chem.MolFromSmiles(s)) if Chem.MolFromSmiles(s) else np.nan)
            elif prop == "LogP":
                df[prop] = df["smiles"].apply(lambda s: Descriptors.MolLogP(Chem.MolFromSmiles(s)) if Chem.MolFromSmiles(s) else np.nan)
            elif prop == "TPSA":
                df[prop] = df["smiles"].apply(lambda s: Descriptors.TPSA(Chem.MolFromSmiles(s)) if Chem.MolFromSmiles(s) else np.nan)
            else:
                continue

            plot_2d_projection(
                umap_result, color=df[prop],
                color_label=prop,
                title=f"UMAP colored by {prop}",
                cmap="viridis",
                output_file=output_dir / f"umap_{prop}.png"
            )

    # ==============================
    #  DBSCAN / Cluster analysis
    # ==============================
    from sklearn.cluster import DBSCAN
    from rdkit import DataStructs
    from rdkit.Chem import Draw
    from rdkit.DataStructs.cDataStructs import CreateFromBitString

    def numpy_to_bitvect(arr: np.ndarray):
        """Convert 0/1 numpy array to RDKit ExplicitBitVect."""
        bitstring = "".join(map(str, arr.astype(int)))
        return CreateFromBitString(bitstring)

    cluster_params = config["dimensionality_reduction_analyses"].get("clustering", {})
    if cluster_params.get("enabled", False):
        logger.info("[DEL] Running DBSCAN clustering on UMAP coordinates...")

        eps = cluster_params.get("eps", 0.5)
        min_samples = cluster_params.get("min_samples", 10)

        db = DBSCAN(eps=eps, min_samples=min_samples)
        cluster_labels = db.fit_predict(umap_result)
        df["umap_cluster"] = cluster_labels

        n_clusters = len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)
        n_noise = list(cluster_labels).count(-1)
        logger.info(f"Found {n_clusters} clusters and {n_noise} noise points")

        cluster_output_dir = output_dir / "clusters"
        cluster_output_dir.mkdir(exist_ok=True)

        # -----------------------------
        # Cluster summary table
        # -----------------------------
        summary_records = []
        for cluster_id, group in df.groupby("umap_cluster"):
            summary_records.append({
                "cluster_id": cluster_id,
                "n_members": len(group),
                "mean_log_enrichment": group["log_enrichment"].mean(),
                "std_log_enrichment": group["log_enrichment"].std()
            })

        cluster_summary_df = (
            pd.DataFrame(summary_records)
            .sort_values("mean_log_enrichment", ascending=False)
        )
        cluster_summary_path = cluster_output_dir / "cluster_summary.csv"
        cluster_summary_df.to_csv(cluster_summary_path, index=False)
        logger.info(f"Saved cluster summary: {cluster_summary_path}")

        # -----------------------------
        # Visualizations
        # -----------------------------
        if cluster_params.get("visualize", False):
            # (1) UMAP colored by cluster
            cmap = plt.cm.get_cmap("tab20", len(set(cluster_labels)))
            plot_2d_projection(
                umap_result,
                color=cluster_labels,
                color_label="Cluster ID",
                title="UMAP colored by DBSCAN cluster",
                cmap=cmap,
                output_file=cluster_output_dir / "umap_clusters.png"
            )

            # (2) Enrichment overlay with cluster IDs
            plt.figure(figsize=(8, 6))
            scatter = plt.scatter(
                umap_result[:, 0], umap_result[:, 1],
                c=df["log_enrichment"], cmap="coolwarm",
                s=40, alpha=0.8, edgecolors="none"
            )
            for clust_id, group in df.groupby("umap_cluster"):
                if clust_id == -1:
                    continue
                centroid = umap_result[group.index].mean(axis=0)
                plt.text(
                    centroid[0], centroid[1], str(clust_id),
                    fontsize=9, fontweight="bold", color="black"
                )
            plt.colorbar(scatter, label="log10(target/control)")
            plt.title("UMAP clusters annotated on enrichment map")
            plt.savefig(
                cluster_output_dir / "umap_cluster_annotations.png",
                dpi=300, bbox_inches="tight"
            )
            plt.close()

            # (3) Enrichment boxplot per cluster
            if cluster_params.get("boxplot", False):
                plt.figure(figsize=(10, 6))
                sns.boxplot(
                    x="umap_cluster",
                    y="log_enrichment",
                    data=df,
                    palette="vlag",
                    hue="umap_cluster",
                    legend=False
                )
                plt.title("Distribution of log10(target/control) per cluster")
                plt.xlabel("Cluster ID")
                plt.ylabel("log10(target/control)")
                plt.xticks(rotation=90)
                plt.tight_layout()
                plt.savefig(
                    cluster_output_dir / "cluster_enrichment_boxplot.png",
                    dpi=300
                )
                plt.close()

        # -----------------------------
        # Representative molecule per cluster
        # -----------------------------
        if cluster_params.get("centroid_images", False):
            centroid_dir = cluster_output_dir / "centroid_molecules"
            centroid_dir.mkdir(exist_ok=True)

            centroid_records = []
            for clust_id, group in df.groupby("umap_cluster"):
                if clust_id == -1 or len(group) < 2:
                    continue

                # Convert stored numpy bit arrays to ExplicitBitVects
                fps = [numpy_to_bitvect(fp) for fp in group["fp_final"]]
                n = len(fps)
                sim_matrix = np.zeros((n, n))
                for i in range(n):
                    for j in range(i + 1, n):
                        sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
                        sim_matrix[i, j] = sim_matrix[j, i] = sim

                mean_sim = sim_matrix.mean(axis=1)
                medoid_idx = np.argmax(mean_sim)
                smiles_central = group.iloc[medoid_idx]["smiles"]
                mol_central = Chem.MolFromSmiles(smiles_central)
                if mol_central is None:
                    continue

                img_path = centroid_dir / f"cluster_{clust_id}_centroid.png"
                img = Draw.MolToImage(mol_central, size=(250, 250))
                img.save(img_path)

                centroid_records.append({
                    "cluster_id": clust_id,
                    "n_members": len(group),
                    "mean_log_enrichment": group["log_enrichment"].mean(),
                    "centroid_smiles": smiles_central,
                    "centroid_image": str(img_path)
                })

            centroid_df = pd.DataFrame(centroid_records)
            centroid_summary_path = cluster_output_dir / "cluster_centroids_summary.csv"
            centroid_df.to_csv(centroid_summary_path, index=False)
            logger.info(f"Saved cluster centroid summary: {centroid_summary_path}")

        # -----------------------------
        # Inter-cluster similarity heatmap
        # -----------------------------
        if cluster_params.get("similarity_heatmap", False):
            centroid_fps = []
            centroid_ids = []
            for clust_id, group in df.groupby("umap_cluster"):
                if clust_id == -1 or len(group) < 2:
                    continue
                fp_central = numpy_to_bitvect(group["fp_final"].iloc[0])
                centroid_fps.append(fp_central)
                centroid_ids.append(clust_id)

            n_centroids = len(centroid_fps)
            if n_centroids > 1:
                sim_mat = np.zeros((n_centroids, n_centroids))
                for i in range(n_centroids):
                    for j in range(i + 1, n_centroids):
                        sim = DataStructs.TanimotoSimilarity(centroid_fps[i], centroid_fps[j])
                        sim_mat[i, j] = sim_mat[j, i] = sim

                plt.figure(figsize=(6, 5))
                sns.heatmap(
                    sim_mat,
                    cmap="viridis",
                    xticklabels=centroid_ids,
                    yticklabels=centroid_ids,
                    cbar_kws={"label": "Tanimoto similarity"}
                )
                plt.title("Inter-cluster Tanimoto similarity")
                plt.tight_layout()
                plt.savefig(
                    cluster_output_dir / "cluster_similarity_heatmap.png",
                    dpi=300
                )
                plt.close()



    # --- Synthon-level enrichment ----------------------------------------------
    if config["dimensionality_reduction_analyses"].get("synthon_enrichment", False):
        for col in ["smiles_a", "smiles_b", "smiles_c"]:
            syn_enr = (
                df.groupby(col)["log_enrichment"]
                .mean()
                .sort_values(ascending=False)
            )
            syn_enr.to_csv(output_dir / f"synthon_enrichment_{col}.csv")
        logger.info("Saved synthon enrichment summaries.")

    # --- Per-synthon analysis ---------------------------------------------------
    if config["dimensionality_reduction_analyses"].get("per_synthon_analysis", False):
        synthon_cols = ["fp_a", "fp_b", "fp_c", "fp_final"]
        for i, synthon_col in enumerate(synthon_cols, start=1):
            logger.info(f"Running per-synthon UMAP: {synthon_col}")
            fp_matrix = np.vstack(df[synthon_col].to_numpy()).astype(float)
            fp_scaled = scaler.fit_transform(fp_matrix)

            umap_model = umap.UMAP(
                n_components=2,
                metric="jaccard",
                n_neighbors=30,
                min_dist=0.1,
                random_state=42 + i  # unique seed per synthon
            )
            umap_syn = umap_model.fit_transform(fp_scaled)

            plot_2d_projection(
                umap_syn,
                color=df["log_enrichment"],
                color_label="log10(target/control)",
                title=f"UMAP ({synthon_col}) colored by enrichment factor",
                output_file=output_dir / f"umap_{synthon_col}.png"
            )

    # --- Save reduced coordinates ----------------------------------------------
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
