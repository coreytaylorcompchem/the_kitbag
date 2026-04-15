import importlib

import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem

import torch
from torch_geometric.loader import DataLoader
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, normalize
from sklearn.decomposition import PCA

from pipeline.task_registry import register_task

from modules.utils.plotting import plot_label_histograms
from modules.utils.splits import scaffold_split

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def quantile_balanced_sampling(
    df,
    label_col="pIC50",
    n_bins=20,
    max_samples_per_bin=None,
    random_state=42,
):
    rng = np.random.default_rng(random_state)

    df = df.copy()
    df["bin"] = pd.qcut(df[label_col], q=n_bins, labels=False, duplicates="drop")

    balanced_idx = []

    for bin_id in df["bin"].unique():
        bin_idx = df[df["bin"] == bin_id].index.to_numpy()

        if max_samples_per_bin:
            n_keep = min(len(bin_idx), max_samples_per_bin)
            chosen = rng.choice(bin_idx, size=n_keep, replace=False)
        else:
            chosen = bin_idx

        balanced_idx.extend(chosen)

    return df.loc[balanced_idx].drop(columns="bin").reset_index(drop=True)


def stratified_balanced_sampling(
    df,
    label_col="pIC50",
    n_bins=20,
    max_samples_per_bin=500,
    random_state=42,
):
    """
    Stratified downsampling to balance dataset by label distribution
    WITHOUT duplication and WITHOUT increasing dataset size.
    """

    rng = np.random.default_rng(random_state)

    # Create bins
    bins = np.linspace(df[label_col].min(), df[label_col].max(), n_bins + 1)
    df = df.copy()
    df["bin"] = pd.cut(df[label_col], bins=bins, include_lowest=True, labels=False)

    balanced_indices = []

    for bin_id in range(n_bins):
        bin_indices = df[df["bin"] == bin_id].index.to_numpy()

        if len(bin_indices) == 0:
            continue

        # Downsample only
        n_keep = min(len(bin_indices), max_samples_per_bin)
        chosen = rng.choice(bin_indices, size=n_keep, replace=False)

        balanced_indices.extend(chosen)

    balanced_df = df.loc[balanced_indices].drop(columns=["bin"]).reset_index(drop=True)

    return balanced_df

@register_task("load_smiles_dataset", category="ADME", description="Load dataset with sampling.")
def load_smiles_dataset(config, context):
    csv_path = config["csv_path"]
    label_col = config.get("label_col", "pIC50")

    df = pd.read_csv(csv_path)
    df_before = df.copy()

    logger.info(f"Dataset size before sampling: {len(df)}")

    sampling_cfg = config.get("sampling")

    if sampling_cfg:
        sample_type = sampling_cfg.get("sample_type")
        n_bins = sampling_cfg.get("n_bins", 20)
        max_samples_per_bin = sampling_cfg.get("max_samples_per_bin", 500)
        random_state = sampling_cfg.get("random_state", 42)

        if sample_type == "balanced":
            logger.info(
                f"Applying stratified balanced sampling: "
                f"{n_bins} bins, max {max_samples_per_bin} per bin"
            )
            df = stratified_balanced_sampling(
                df,
                label_col=label_col,
                n_bins=n_bins,
                max_samples_per_bin=max_samples_per_bin,
                random_state=random_state,
            )

        elif sample_type == "quantile":
            logger.info(
                f"Applying quantile balanced sampling: "
                f"{n_bins} bins, max {max_samples_per_bin} per bin"
            )
            df = quantile_balanced_sampling(
                df,
                label_col=label_col,
                n_bins=n_bins,
                max_samples_per_bin=max_samples_per_bin,
                random_state=random_state,
            )

        else:
            logger.info("Sampling block present but no valid sample_type specified")

        # Optional diagnostic plot
        if sampling_cfg.get("plot_distributions", False):
            plot_dir = sampling_cfg.get("plot_dir", "outputs/cyp3a4/sampling")
            plot_path = plot_label_histograms(
                df_before,
                df,
                label_col=label_col,
                out_dir=plot_dir,
                title_suffix=f" ({sample_type})",
            )
            logger.info(f"Saved label distribution plot to {plot_path}")

    else:
        logger.info("No sampling configuration found; using full dataset")

    context["dataframe"] = df
    logger.info(f"Dataset size after sampling: {len(df)}")

    return {"dataframe": df}

def ic50_to_pic50(x):
    return -np.log10(x * 1e-9)  # assume values are in nM

@register_task("transform_labels", category="ADME")
def transform_labels(config, context):

    label_cols = config["label_cols"] # from the yaml

    df = context["dataframe"]
    task_cols = config["label_cols"]

    for col in task_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')

        if "IC50" in col:
            df[col] = df[col].apply(
                lambda x: ic50_to_pic50(x) if pd.notna(x) and x > 0 else np.nan
            )

        elif col in ["LogP", "LogD"]:
            # realistic ADME bounds
            lower, upper = -5, 10

            outlier_mask = (df[col] < lower) | (df[col] > upper)

            n_outliers = outlier_mask.sum()

            if n_outliers > 0:
                logger.warning(
                    f"{col}: removing {n_outliers} outliers outside [{lower}, {upper}]"
                )
            
            df.loc[outlier_mask, col] = np.nan # removes the outliers

        elif col == "Solubility":
            # don't log transform LogS
            # TODO: automatically detect entirely negative data
            pass
        
        elif "Microsomal Stability" in col:
            # Clint-specific handling

            # Remove invalid values
            valid_mask = df[col].notna()

            invalid_mask = valid_mask & (
                (df[col] <= 0) | (~np.isfinite(df[col]))
            )

            n_invalid = invalid_mask.sum()
            if n_invalid > 0:
                logger.warning(f"{col}: removing {n_invalid} invalid (<=0 or inf) values")

            df.loc[invalid_mask, col] = np.nan

            upper = df[col].quantile(0.99)
            df.loc[df[col] > upper, col] = upper

            # Log transform 
            df.loc[df[col].notna(), col] = np.log(df.loc[df[col].notna(), col])

        else:
            df[col] = np.log1p(df[col])

    scalers = {}

    for col in task_cols:
        mask = df[col].notna()

        if mask.sum() == 0:
            continue

        scaler = StandardScaler()

        scaled = scaler.fit_transform(
            df.loc[mask, col].values.reshape(-1, 1)
        ).flatten()

        df.loc[mask, col] = scaled

        scalers[col] = scaler

    # Save scalers for later use
    context["label_scalers"] = scalers

    context["task_names"] = label_cols
    context["dataframe"] = df

    return context

# @register_task("compute_mtl_features", category="ADME", description="Compute global features and fingerprints for MTL")
# def compute_mtl_features(config, context):
#     df = context["dataframe"]
#     smiles_list = df["smiles"].tolist()

#     # =============================
#     # 1. GLOBAL FEATURES
#     # =============================
#     from models.adme.mtl_adme.featurisation import process_single_mol

#     inputs = list(enumerate(smiles_list))

#     global_feats = Parallel(n_jobs=config.get("n_jobs", -1))(
#         delayed(process_single_mol)((i, smi))
#         for i, smi in tqdm(inputs, desc="Global features")
#     )

#     global_feats = np.array(global_feats, dtype=np.float32)

#     scaler = StandardScaler()
#     scaled_global_feats = scaler.fit_transform(global_feats)

#     # =============================
#     # 2. FINGERPRINTS + PCA
#     # =============================
#     from models.adme.mtl_adme.featurisation import compute_morgan_fp

#     all_fps = []

#     for smi in tqdm(smiles_list, desc="Fingerprints"):
#         mol = Chem.MolFromSmiles(smi)
#         if mol is None:
#             all_fps.append(np.zeros(1024, dtype=np.float32))
#             continue

#         fp_ecfp = compute_morgan_fp(mol, n_bits=512)
#         fp_torsion = np.array(
#             AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect(mol, nBits=512)
#         )

#         fp = np.concatenate([fp_ecfp, fp_torsion])
#         all_fps.append(fp)

#     all_fps = np.array(all_fps, dtype=np.float32)

#     # PCA → 512 dims (same as notebook)
#     pca = PCA(n_components=512)
#     fps_pca = pca.fit_transform(all_fps)

#     # Normalize
#     fps_pca = normalize(fps_pca, norm="l2", axis=1)

#     # =============================
#     # Save to context
#     # =============================
#     context["scaled_global_feats"] = scaled_global_feats
#     context["fps_pca"] = fps_pca

#     return {
#         "scaled_global_feats": scaled_global_feats,
#         "fps_pca": fps_pca
#     }

@register_task("featurise_smiles", category="ADME", description="Featurise SMILES for graph-based models.")
def featurise_smiles(config, context):
    df = context["dataframe"]

    smiles_col = config.get("smiles_col", "smiles")

    # Retrieve data and label cols for multi-task
    label_cols = config.get("label_cols")
    label_col = config.get("label_col")

    if label_cols is None and label_col is None:
        raise ValueError("Provide either 'label_col' or 'label_cols'")

    # Load featuriser dynamically
    feat_cfg = config["featuriser"]
    module = importlib.import_module(feat_cfg["module"])
    featuriser = getattr(module, feat_cfg["function"])

    if hasattr(module, "prepare_features"):
        module.prepare_features(df, smiles_col=smiles_col)

    graphs = []

    # MULTI-TASK
    if label_cols is not None:
        labels = df[label_cols].values.astype(np.float32)

        for i, (smi, y_vec) in enumerate(
            tqdm(zip(df[smiles_col], labels),
                total=len(df),
                desc="Featurising molecules")
        ):
            g = featuriser(smi, label=y_vec, idx=i)
            if g is not None:
                graphs.append(g)

        num_tasks = len(label_cols)

    # SINGLE-TASK (backward compatible)
    else:
        for i, (smi, y) in enumerate(
            tqdm(zip(df[smiles_col], df[label_col]),
                total=len(df),
                desc="Featurising molecules")
        ):
            g = featuriser(smi, label=y, idx=i)
            if g is not None:
                graphs.append(g)

        num_tasks = 1

    if not graphs:
        raise ValueError("Featurisation produced zero graphs.")

    # Infer dimensions
    first_graph = graphs[0]

    global_feat = getattr(first_graph, "global_features", None)
    fp = getattr(first_graph, "fp", None)

    context.update({
        "graphs": graphs,
        "input_dim": first_graph.x.shape[1],
        "edge_dim": first_graph.edge_attr.shape[1],
        "global_feat_dim": global_feat.shape[-1] if global_feat is not None else 0,
        "fp_dim": fp.shape[-1] if fp is not None else 0,
        "num_tasks": num_tasks
    })

    return context

@register_task("split_data", category="ADME", description="Perform train/test/val splits.")
def split_data(config, context):
    graphs = context["graphs"]
    smiles = context["dataframe"]["smiles"].tolist()

    test_size = config.get("val_size", 0.2)
    batch_size = config.get("batch_size", 32)

    train_list, val_list = scaffold_split(
        graphs,
        smiles,
        val_fraction=config.get("val_size", 0.2),
        seed=42
    )

    train_loader = DataLoader(train_list, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_list, batch_size=batch_size, shuffle=False)

    context.update({
        "train_loader": train_loader,
        "val_loader": val_loader
    })

    return {
        "train_loader": train_loader,
        "val_loader": val_loader
    }

@register_task("load_adme_model_spec", category="ADME", description="Instantiate Pytorch ADME models.")
def load_adme_model_spec(config, context):
    if "module" not in config or "class" not in config:
        raise KeyError(
            "load_adme_model_spec requires 'module' and 'class' in YAML config. "
            "Check that the task name matches the YAML block."
        )
    
    module_path = config["module"]
    class_name = config["class"]

    module = importlib.import_module(module_path)
    ModelClass = getattr(module, class_name)

    context["model_class"] = ModelClass
    return {"model_class": ModelClass}

@register_task("train_adme_model", category="ADME", description="Train Pytorch model.")
def train_adme_model(config, context):
    trainer_cfg = config["trainer"]
    module = importlib.import_module(trainer_cfg["module"])
    train_fn = getattr(module, trainer_cfg["function"])

    context["device"] = torch.device(config.get("device", "cuda" if torch.cuda.is_available() else "cpu"))

    # Call trainer
    result = train_fn(context=context, config=config)

    # Handle either dict or direct model
    if isinstance(result, dict):
        trained_model = result.get("model")
    else:
        trained_model = result

    context.update({
        "trained_model": trained_model,
        "training_summary": result,
    })

    return context

@register_task("evaluate_adme_model", category="ADME", description="Evaluate trained model.")
def evaluate_model(config, context):
    """
    Generic evaluation task.
    Calls a model-specific evaluation function defined in trainer/evaluation module.
    """
    eval_cfg = config.get("evaluator")
    if eval_cfg is None:
        # Fallback: try using module/function from config directly
        module_path = config.get("module", "models.adme.cyp3a4.evaluation")
        function_name = config.get("function", "evaluate")
    else:
        module_path = eval_cfg["module"]
        function_name = eval_cfg["function"]

    # Import evaluation function dynamically
    module = importlib.import_module(module_path)
    eval_fn = getattr(module, function_name)

    # Call the evaluation function with context and config
    result = eval_fn(context=context, config=config)

    # Update context with any results returned by the evaluation function
    if isinstance(result, dict):
        context.update(result)

    return context


