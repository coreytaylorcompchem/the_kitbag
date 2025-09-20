import torch
import pandas as pd
from pathlib import Path
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

# Import featurization and model classes
from models.featurisation_cyp import mol_to_graph as mol_to_graph_cyp
from models.featurisation_logd import mol_to_graph as mol_to_graph_logd  # assuming this exists
from models.model_defs import GINRegressor, GATv2Regressor

import torch_geometric
from torch_geometric.data import Data

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def get_available_adme_models():
    model_dir = Path(__file__).resolve().parent.parent / "models"
    return {
        "cyp3a4": {
            "path": model_dir / "cyp_3a4_gin.pt",
            "class": GINRegressor,
            "featuriser": mol_to_graph_cyp
        },
        "logd": {
            "path": model_dir / "logd_gat.pt",
            "class": GATv2Regressor,
            "featuriser": mol_to_graph_logd
        }
    }

@register_task("adme_prediction", description="Predict ADME properties (hERG, LogD, CYP3A4) from SMILES using PyTorch models")
def adme_prediction(config, data=None):
    requested_models = config.get("adme_models", [])
    if not isinstance(requested_models, list):
        raise ValueError("'adme_models' should be a list of model names.")

    all_models = get_available_adme_models()

    if not requested_models:
        model_configs = all_models
    else:
        model_configs = {}
        for name in requested_models:
            if name not in all_models:
                raise ValueError(f"Unknown ADME model: '{name}'")
            model_configs[name] = all_models[name]

    input_file = config.get("input_file")
    if not input_file:
        raise ValueError("Missing 'input_file' in config.")

    df = pd.read_csv(input_file)

    if "smiles" not in df.columns:
        raise ValueError("Input file must contain a 'smiles' column.")

    df["smiles"] = df["smiles"].astype(str).str.strip()
    df = df[df["smiles"].notna() & df["smiles"].str.len() > 0].copy()

    result_df = df[["smiles"]].copy()

    device = torch.device("cpu" if torch.cuda.is_available() else "cpu")
    logger.debug(f"Using device: {device}")

    for model_name, info in model_configs.items():
        model_path = info["path"]
        model_class = info["class"]
        logger.debug(f"Loading model '{model_name}' from {model_path}")

        try:
            checkpoint = torch.load(model_path, map_location=device)
            model_args = {
                "input_dim": checkpoint["input_dim"],
                "hidden_dim": checkpoint.get("hidden_dim", 128),
                "global_feat_dim": checkpoint.get("global_feat_dim", 0)
            }
            if "edge_dim" in checkpoint and "edge_dim" in model_class.__init__.__code__.co_varnames:
                model_args["edge_dim"] = checkpoint["edge_dim"]

            model = model_class(**model_args)
            model.load_state_dict(checkpoint["model_state_dict"])
            model.to(device)
            model.eval()
        except Exception as e:
            logger.error(f"Failed to load model '{model_name}': {e}")
            result_df[model_name] = None
            continue

        preds = []
        try:
            preds = []
            featuriser = info.get("featuriser")  # use model-specific featuriser
            for smi in df["smiles"]:
                graph_data = featuriser(smi)
                if graph_data is None:
                    logger.error(f"[model] Prediction failed for SMILES: '{smi}' — featuriser returned None")
                    preds.append(None)
                    continue
                
                logger.debug(f"graph_data type: {type(graph_data)} for SMILES: {smi}")
                graph_data = graph_data.to(device)

                with torch.no_grad():
                    output = model(graph_data)
                    preds.append(output.item() if isinstance(output, torch.Tensor) else float(output))

            result_df[model_name] = preds

        except Exception as e:
            logger.error(f"Prediction failed for model '{model_name}': {e}")
            result_df[model_name] = None

    return (Path(input_file).stem, result_df)
