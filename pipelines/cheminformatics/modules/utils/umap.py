import numpy as np
import warnings

import umap
import hdbscan
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

def add_umap_clusters(
    df,
    feature_cols,
    n_neighbors=30,
    min_dist=0.1,
    n_components=2,
    cluster_method="hdbscan",
    random_state=42,
    min_points_for_umap=10,
):
    """
    Adds UMAP embedding + cluster labels for plotting only.
    Safe against small-N and degenerate cases.
    """

    # ----------------------------
    # Feature prep
    # ----------------------------
    X = df[feature_cols].replace([np.inf, -np.inf], np.nan).dropna()

    if len(X) < min_points_for_umap:
        # Too few points → skip UMAP entirely
        df["umap_x"] = np.nan
        df["umap_y"] = np.nan
        df["cluster_id"] = -1

        warnings.warn(
            f"[UMAP] Skipped: only {len(X)} points (min={min_points_for_umap})."
        )
        return df

    valid_idx = X.index
    X = X.values

    X_scaled = StandardScaler().fit_transform(X)

    # ----------------------------
    # Adaptive neighbors
    # ----------------------------
    n_neighbors_eff = min(n_neighbors, len(X_scaled) - 1)
    if n_neighbors_eff < 2:
        df["umap_x"] = np.nan
        df["umap_y"] = np.nan
        df["cluster_id"] = -1
        return df

    try:
        reducer = umap.UMAP(
            n_neighbors=n_neighbors_eff,
            min_dist=min_dist,
            n_components=n_components,
            random_state=random_state,
        )

        embedding = reducer.fit_transform(X_scaled)

    except Exception as e:
        df["umap_x"] = np.nan
        df["umap_y"] = np.nan
        df["cluster_id"] = -1

        warnings.warn(f"[UMAP] Failed, falling back. Reason: {e}")
        return df

    df.loc[valid_idx, "umap_x"] = embedding[:, 0]
    df.loc[valid_idx, "umap_y"] = embedding[:, 1]

    # ----------------------------
    # Clustering
    # ----------------------------
    try:
        if cluster_method == "hdbscan":
            clusterer = hdbscan.HDBSCAN(
                min_cluster_size=max(5, len(embedding) // 20),
                min_samples=3
            )
            labels = clusterer.fit_predict(embedding)

        elif cluster_method == "kmeans":
            
            k = max(2, min(10, len(embedding) // 5))
            labels = KMeans(n_clusters=k, random_state=random_state).fit_predict(embedding)

        else:
            raise ValueError(f"Unknown cluster_method: {cluster_method}")

        df.loc[valid_idx, "cluster_id"] = labels

    except Exception as e:
        df["cluster_id"] = -1
        warnings.warn(f"[Clustering] Failed, assigning noise. Reason: {e}")

    return df

