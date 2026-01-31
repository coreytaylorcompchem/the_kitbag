import importlib

import pandas as pd
import numpy as np

import torch
from torch_geometric.loader import DataLoader
from sklearn.model_selection import train_test_split

from pipeline.task_registry import register_task
from modules.utils.plotting import plot_label_histograms

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

# def stratified_balanced_sampling(df, label_col='pIC50', n_bins=20, max_samples_per_bin=500, random_state=42):
#     """
#     Stratified sampling to balance dataset by label distribution.

#     Args:
#         df (pd.DataFrame): DataFrame containing labelled data.
#         label_col (str): Name of the label column (e.g., 'pIC50').
#         n_bins (int): Number of bins to discretize the label.
#         max_samples_per_bin (int): Maximum number of samples to keep per bin.
#         random_state (int): Random seed for reproducibility.

#     Returns:
#         pd.DataFrame: Balanced DataFrame sampled from input.
#     """
#     np.random.seed(random_state)
    
#     bins = np.linspace(df[label_col].min(), df[label_col].max(), n_bins + 1)
#     df['bin'] = pd.cut(df[label_col], bins=bins, include_lowest=True, labels=False)

#     balanced_indices = []

#     for bin_id in range(n_bins):
#         bin_indices = df[df['bin'] == bin_id].index.values
#         n_samples = len(bin_indices)

#         if n_samples == 0:
#             continue
        
#         if n_samples > max_samples_per_bin:
#             chosen = np.random.choice(bin_indices, max_samples_per_bin, replace=False)
#         else:
#             chosen = np.random.choice(bin_indices, max_samples_per_bin, replace=True)

#         balanced_indices.extend(chosen)

#     balanced_df = df.loc[balanced_indices].reset_index(drop=True)

#     balanced_df = balanced_df.drop(columns=['bin'])

#     return balanced_df

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


@register_task("featurise_smiles", category="ADME", description="Featurise SMILES for graph-based models.")
def featurise_smiles(config, context):
    df = context["dataframe"]

    smiles_col = config.get("smiles_col", "smiles")
    label_col = config.get("label_col", "pIC50")

    # Load featuriser dynamically
    feat_cfg = config["featuriser"]
    module = importlib.import_module(feat_cfg["module"])
    featuriser = getattr(module, feat_cfg["function"])

    graphs = []
    for smi, y in zip(df[smiles_col], df[label_col]):
        g = featuriser(smi, label=y)
        if g is not None:
            graphs.append(g)

    if not graphs:
        raise ValueError("Featurisation produced zero graphs.")

    # Infer dimensions for downstream tasks
    context.update({
        "graphs": graphs,
        "input_dim": graphs[0].x.shape[1],
        "global_feat_dim": getattr(graphs[0], "global_features", None).shape[0]
    })

    return context

@register_task("split_data", category="ADME", description="Perform train/test/val splits.")
def split_data(config, context):
    graphs = context["graphs"]

    test_size = config.get("val_size", 0.2)
    batch_size = config.get("batch_size", 32)

    train_list, val_list = train_test_split(
        graphs,
        test_size=test_size,
        random_state=config.get("seed", 42)
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


