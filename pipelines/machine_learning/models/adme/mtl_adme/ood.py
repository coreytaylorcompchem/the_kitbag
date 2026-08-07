import numpy as np
import pandas as pd

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors
from sklearn.preprocessing import RobustScaler

from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Cluster import Butina

DESCRIPTOR_FUNCS = {
    "MolWt": Descriptors.MolWt,
    "MolLogP": Descriptors.MolLogP,
    "TPSA": Descriptors.TPSA,
    "NumHDonors": Descriptors.NumHDonors,
    "NumHAcceptors": Descriptors.NumHAcceptors,
    "NumRotatableBonds": Descriptors.NumRotatableBonds,
    "FractionCSP3": Descriptors.FractionCSP3,
    "RingCount": Descriptors.RingCount,
    "HeavyAtomCount": Descriptors.HeavyAtomCount,
}

def murcko_scaffold(smiles):
    mol = mol_from_smiles(smiles)

    if mol is None:
        return None

    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaffold)
    except Exception:
        return None

def compute_scaffold_features(query_smiles, train_smiles):

    train_scaffolds = {
        murcko_scaffold(s)
        for s in train_smiles
    }

    rows = []

    for smi in query_smiles:
        scaffold = murcko_scaffold(smi)
        rows.append({
            "murcko_scaffold": scaffold,
            "known_scaffold":
                scaffold in train_scaffolds
                if scaffold is not None
                else False
        })

    return pd.DataFrame(rows)

def mol_from_smiles(smiles):
    if smiles is None:
        return None
    try:
        mol = Chem.MolFromSmiles(str(smiles))
    except Exception:
        mol = None

    return mol


def morgan_fp_from_smiles(smiles, radius=2, n_bits=2048):
    mol = mol_from_smiles(smiles)

    if mol is None:
        return None

    return AllChem.GetMorganFingerprintAsBitVect(
        mol,
        radius,
        nBits=n_bits
    )

def fingerprints_from_smiles(
    smiles_list,
    radius=2,
    n_bits=2048,
):
    """
    Generate Morgan fingerprints while preserving indexing.
    """

    fps = []

    for smi in smiles_list:
        fps.append(
            morgan_fp_from_smiles(
                smi,
                radius=radius,
                n_bits=n_bits,
            )
        )

    return fps

def cluster_training_set(
    train_fps,
    cutoff=0.4,
):
    """
    Cluster the training fingerprints using Butina.

    Returns
    -------
    clusters : tuple(tuple(int))
        Cluster membership.

    cluster_sizes : ndarray
        Size of each cluster.

    molecule_cluster : ndarray
        Cluster index for every molecule.
    """

    valid = [
        fp
        for fp in train_fps
        if fp is not None
    ]

    if len(valid) == 0:
        raise ValueError("No valid fingerprints.")

    dists = []

    for i in range(1, len(valid)):
        sims = DataStructs.BulkTanimotoSimilarity(
            valid[i],
            valid[:i],
        )

        dists.extend(
            [1.0 - s for s in sims]
        )

    clusters = Butina.ClusterData(
        dists,
        len(valid),
        cutoff,
        isDistData=True,
    )

    cluster_sizes = np.zeros(
        len(clusters),
        dtype=int,
    )

    molecule_cluster = np.zeros(
        len(valid),
        dtype=int,
    )

    for cid, members in enumerate(clusters):

        cluster_sizes[cid] = len(members)

        for m in members:
            molecule_cluster[m] = cid

    return (
        clusters,
        cluster_sizes,
        molecule_cluster,
    )

def assign_query_cluster(
    query_fp,
    train_fps,
    molecule_cluster,
):
    """
    Assign a query molecule to the cluster containing its nearest
    neighbour.
    """

    if query_fp is None:
        return None

    sims = DataStructs.BulkTanimotoSimilarity(
        query_fp,
        train_fps,
    )

    idx = int(np.argmax(sims))

    return molecule_cluster[idx]

def compute_cluster_features(
    query_smiles,
    train_smiles,
    cutoff=0.4,
):
    """
    Compute cluster-based applicability-domain features.
    """

    train_fps = [
        fp
        for fp in fingerprints_from_smiles(train_smiles)
        if fp is not None
    ]

    (
        clusters,
        cluster_sizes,
        molecule_cluster,
    ) = cluster_training_set(
        train_fps,
        cutoff=cutoff,
    )

    rows = []

    for smi in query_smiles:
        fp = morgan_fp_from_smiles(smi)
        cid = assign_query_cluster(
            fp,
            train_fps,
            molecule_cluster,
        )

        if cid is None:
            rows.append({
                "cluster_id": np.nan,
                "cluster_size": np.nan,
            })

            continue

        rows.append({
            "cluster_id": cid,
            "cluster_size": cluster_sizes[cid],
        })

    return pd.DataFrame(rows)


def compute_nn_similarity_features(
    query_smiles,
    train_smiles,
    radius=2,
    n_bits=2048,
    top_k=5,
    similarity_thresholds=(0.4, 0.5, 0.6, 0.7),
):
    """
    Compute nearest-neighbour and local-density similarity features.

    Memory behaviour:
      - Stores train fingerprints once.
      - Processes one query molecule at a time.
      - Does NOT build a query x train similarity matrix.
    """

    train_fps = [
        morgan_fp_from_smiles(
            smi,
            radius=radius,
            n_bits=n_bits,
        )
        for smi in train_smiles
    ]

    valid_train_fps = [
        fp for fp in train_fps
        if fp is not None
    ]

    if len(valid_train_fps) == 0:
        raise ValueError(
            "No valid training fingerprints were generated. "
            "Check train_smiles and RDKit parsing."
        )

    rows = []
    n_train = len(valid_train_fps)

    for smi in query_smiles:
        qfp = morgan_fp_from_smiles(
            smi,
            radius=radius,
            n_bits=n_bits,
        )

        if qfp is None:
            row = {
                "nearest_train_tanimoto": np.nan,
                f"mean_top{top_k}_train_tanimoto": np.nan,
            }

            for thr in similarity_thresholds:
                row[f"n_train_neighbors_ge_{thr}"] = np.nan
                row[f"fraction_train_neighbors_ge_{thr}"] = np.nan

            rows.append(row)
            continue

        sims = np.asarray(
            DataStructs.BulkTanimotoSimilarity(
                qfp,
                valid_train_fps,
            ),
            dtype=np.float32,
        )

        if sims.size == 0:
            row = {
                "nearest_train_tanimoto": np.nan,
                f"mean_top{top_k}_train_tanimoto": np.nan,
            }

            for thr in similarity_thresholds:
                row[f"n_train_neighbors_ge_{thr}"] = np.nan
                row[f"fraction_train_neighbors_ge_{thr}"] = np.nan

            rows.append(row)
            continue

        k_eff = min(top_k, sims.size)

        if k_eff == sims.size:
            topk = np.sort(sims)[::-1]
        else:
            topk = np.partition(
                sims,
                -k_eff,
            )[-k_eff:]

        nearest = float(np.max(topk))
        mean_topk = float(np.mean(topk))

        row = {
            "nearest_train_tanimoto": nearest,
            f"mean_top{top_k}_train_tanimoto": mean_topk,
        }

        for thr in similarity_thresholds:
            n_ge = int(np.sum(sims >= float(thr)))

            row[f"n_train_neighbors_ge_{thr}"] = n_ge
            row[f"fraction_train_neighbors_ge_{thr}"] = n_ge / n_train

        rows.append(row)

    return pd.DataFrame(rows)

def compute_physchem_descriptor_frame(smiles_list, descriptors):
    """
    Compute RDKit physicochemical descriptors for a list of SMILES.

    Invalid molecules get NaN values for all requested descriptors.
    """

    rows = []

    for smi in smiles_list:
        mol = mol_from_smiles(smi)

        if mol is None:
            rows.append({
                d: np.nan
                for d in descriptors
            })
            continue

        row = {}

        for d in descriptors:
            if d not in DESCRIPTOR_FUNCS:
                raise ValueError(
                    f"Unknown descriptor '{d}'. "
                    f"Available descriptors: {sorted(DESCRIPTOR_FUNCS)}"
                )

            try:
                row[d] = DESCRIPTOR_FUNCS[d](mol)
            except Exception:
                row[d] = np.nan

        rows.append(row)

    return pd.DataFrame(rows)


def compute_physchem_distance_features(
    query_smiles,
    train_smiles,
    descriptors,
):
    """
    Compute robust-scaled Euclidean distance from the training-set
    physicochemical centre.
    """

    train_desc = compute_physchem_descriptor_frame(
        train_smiles,
        descriptors
    )

    query_desc = compute_physchem_descriptor_frame(
        query_smiles,
        descriptors
    )

    train_medians = train_desc.median(numeric_only=True)

    train_desc = train_desc.fillna(train_medians)
    query_desc = query_desc.fillna(train_medians)

    scaler = RobustScaler()

    train_scaled = scaler.fit_transform(train_desc.values)
    query_scaled = scaler.transform(query_desc.values)

    train_center = np.median(train_scaled, axis=0)

    distances = np.sqrt(
        np.sum(
            (query_scaled - train_center) ** 2,
            axis=1
        )
    )

    return (
        pd.DataFrame({
            "physchem_robust_distance": distances
        }),
        train_scaled,
    )


def minmax_scale_series(values):
    values = np.asarray(values, dtype=float)

    out = np.full_like(
        values,
        np.nan,
        dtype=float
    )

    mask = np.isfinite(values)

    if mask.sum() == 0:
        return out

    vmin = np.nanmin(values[mask])
    vmax = np.nanmax(values[mask])

    if np.isclose(vmin, vmax):
        out[mask] = 0.0
    else:
        out[mask] = (
            values[mask] - vmin
        ) / (
            vmax - vmin
        )

    return out

def percentile_scale_series(values, reference_values):
    """
    Scale values relative to a reference distribution.

    Returns percentile rank (0-1).

    0 = similar to reference centre
    1 = extreme relative to reference
    """

    values = np.asarray(values, dtype=float)
    reference_values = np.asarray(reference_values, dtype=float)

    out = np.full_like(
        values,
        np.nan,
        dtype=float
    )

    reference_values = reference_values[
        np.isfinite(reference_values)
    ]

    if len(reference_values) == 0:
        return out

    for i, v in enumerate(values):

        if not np.isfinite(v):
            continue

        out[i] = (
            np.sum(reference_values <= v)
            / len(reference_values)
        )

    return out

def compute_composite_ood_score(
    feature_df,
    density_col="n_train_neighbors_ge_0.6",
):
    """
    Combine multiple OOD indicators into one score.

    Higher = further outside training domain.

    Components:
      - low nearest-neighbour similarity
      - high physicochemical distance
      - low local neighbour density
      - novel Murcko scaffold
    """

    similarity = feature_df["nearest_train_tanimoto"].values
    physchem = feature_df["physchem_robust_distance"].values

    if density_col in feature_df.columns:
        density = feature_df[density_col].values
    else:
        density = np.full(len(feature_df), np.nan)

    # High similarity is good, so invert it.
    similarity_component = 1.0 - minmax_scale_series(similarity)

    # Large physchem distance is bad.
    physchem_component = minmax_scale_series(physchem)

    # High neighbour density is good, so invert it.
    density_component = 1.0 - minmax_scale_series(density)

    components = [
        similarity_component,
        physchem_component,
        density_component,
    ]

    if "known_scaffold" in feature_df.columns:
        scaffold_component = np.where(
            feature_df["known_scaffold"].astype(bool).values,
            0.0,
            1.0,
        )
        components.append(scaffold_component)

    component_matrix = np.vstack(components)

    return np.nanmean(
        component_matrix,
        axis=0,
    )

def assign_ood_bins(scores):
    """
    Convert continuous OOD scores into interpretable categories.
    """

    bins = []

    for score in scores:

        if np.isnan(score):
            bins.append("unknown")

        elif score < 0.33:
            bins.append("in_domain")

        elif score < 0.66:
            bins.append("moderate_ood")

        else:
            bins.append("high_ood")

    return bins

def compute_ood_features(
    query_smiles,
    train_smiles,
    config,
):
    """
    Compute applicability-domain features for every query molecule.

    This version is suitable for large datasets because it does not use
    all-pairs Butina clustering.
    """

    similarity_thresholds = config.get(
        "similarity_thresholds",
        (0.4, 0.5, 0.6, 0.7),
    )

    similarity_df = compute_nn_similarity_features(
        query_smiles=query_smiles,
        train_smiles=train_smiles,
        radius=config.get("radius", 2),
        n_bits=config.get("n_bits", 2048),
        top_k=config.get("top_k", 5),
        similarity_thresholds=similarity_thresholds,
    )

    scaffold_df = compute_scaffold_features(
        query_smiles=query_smiles,
        train_smiles=train_smiles,
    )

    physchem_df, _ = compute_physchem_distance_features(
        query_smiles=query_smiles,
        train_smiles=train_smiles,
        descriptors=config.get(
            "physchem_descriptors",
            [
                "MolWt",
                "MolLogP",
                "TPSA",
                "NumHDonors",
                "NumHAcceptors",
                "NumRotatableBonds",
                "FractionCSP3",
            ],
        ),
    )

    feature_df = pd.concat(
        [
            similarity_df,
            scaffold_df,
            physchem_df,
        ],
        axis=1,
    )

    density_threshold = config.get(
        "density_threshold",
        0.6,
    )

    density_col = f"n_train_neighbors_ge_{density_threshold}"

    feature_df["ood_score"] = compute_composite_ood_score(
        feature_df,
        density_col=density_col,
    )

    feature_df["ood_percentile"] = (
        pd.Series(feature_df["ood_score"])
        .rank(pct=True)
        .values
    )

    feature_df["ood_bin"] = assign_ood_bins(
        feature_df["ood_score"].values
    )

    return feature_df
