import torch
import pandas as pd
from pathlib import Path
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

# Import featurization and model classes
from models.featurisation_cyp import mol_to_graph
from models.model_defs import GINRegressor
# TODO: Import other model classes as needed (e.g., HERGClassifier, LogDModel)

logger = setup_logger(__name__, debug_mode=False, simple_format=True)


def get_available_adme_models():
    """Map model names to file paths and architecture classes."""
    model_dir = Path(__file__).resolve().parent.parent / "models"
    return {
        "cyp3a4": {
            "path": model_dir / "cyp_3a4_gin.pt",
            "class": GINRegressor
        },
        # Add others here as needed:
        # "herg": {"path": model_dir / "herg_gnn.pt", "class": HERGClassifier},
        # "logd": {"path": model_dir / "logd_gat.pt", "class": LogDModel},
    }


@register_task("adme_prediction", description="Predict ADME properties (hERG, LogD, CYP3A4) from SMILES using PyTorch models")
def adme_prediction(config, data=None):
    requested_models = config.get("adme_models", [])
    if not isinstance(requested_models, list):
        raise ValueError("'adme_models' should be a list of model names.")

    all_models = get_available_adme_models()

    if not requested_models:
        # Default to all available models
        model_configs = all_models
    else:
        model_configs = {}
        for name in requested_models:
            if name not in all_models:
                raise ValueError(f"Unknown ADME model: '{name}'")
            model_configs[name] = all_models[name]

    # Load input directly from yaml.

    input_file = config.get("input_file")
    if not input_file:
        raise ValueError("Missing 'input_file' in config.")

    df = pd.read_csv(input_file)

    ### Not quite working and not sure why, as it dumps columns when making ADME predictions. Also it's slow.
    # To investigate.

    # input_file = config.get("input_file", None)

    # if isinstance(data, dict):
    #     first_key = next(iter(data))
    #     df = data[first_key]
    #     if not isinstance(df, pd.DataFrame):
    #         raise ValueError(f"Expected DataFrame for key {first_key}, got {type(df)}")
    # else:
    #     if not input_file:
    #         raise ValueError("Missing 'input_file' in config")
    #     df = pd.read_csv(input_file)

    # if df.empty:
    #     logger.info("Input DataFrame is empty after filtering — skipping ADME prediction")
    #     # Return empty DataFrame with smiles and model columns
    #     model_names = config.get("adme_models", [])
    #     if isinstance(model_names, list):
    #         model_columns = model_names
    #     elif isinstance(model_names, dict):
    #         model_columns = list(model_names.keys())
    #     else:
    #         model_columns = []

    #     empty_result_df = pd.DataFrame(columns=["smiles"] + model_columns)
    #     # Return safe stem name or fallback string
    #     stem_name = Path(input_file).stem if input_file else "chunk_empty"
    #     return (stem_name, empty_result_df)

    if "smiles" not in df.columns:
        raise ValueError("Input file must contain a 'smiles' column.")

    df["smiles"] = df["smiles"].astype(str).str.strip()
    df = df[df["smiles"].notna() & df["smiles"].str.len() > 0].copy()

    result_df = df[["smiles"]].copy()
    device = torch.device("cpu" if torch.cuda.is_available() else "cpu")

    for model_name, info in model_configs.items():
        model_path = info["path"]
        model_class = info["class"]
        logger.debug(f"Loading model '{model_name}' from {model_path}")

        try:
            checkpoint = torch.load(model_path, map_location=device)
            model = model_class(
                input_dim=checkpoint["input_dim"],
                hidden_dim=checkpoint["hidden_dim"],
                global_feat_dim=checkpoint.get("global_feat_dim", 0)
            )
            model.load_state_dict(checkpoint["model_state_dict"])
            model.to(device)
            model.eval()
        except Exception as e:
            logger.error(f"Failed to load model '{model_name}': {e}")
            result_df[model_name] = None
            continue

        try:
            preds = []
            for smi in df["smiles"]:
                graph_data = mol_to_graph(smi)
                if graph_data is None:
                    preds.append(None)
                    continue

                graph_data = graph_data.to(device)
                with torch.no_grad():
                    output = model(graph_data)
                    preds.append(output.item() if isinstance(output, torch.Tensor) else float(output))

            result_df[model_name] = preds

        except Exception as e:
            logger.error(f"Prediction failed for model '{model_name}': {e}")
            result_df[model_name] = None

    return (Path(input_file).stem, result_df)