import torch
import pandas as pd
from pathlib import Path
from torch_geometric.loader import DataLoader
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

# Import featurization and model classes
from models.featurisation_cyp import mol_to_graph as mol_to_graph_cyp
from models.featurisation_logd import mol_to_graph as mol_to_graph_logd
from models.featurisation_herg import mol_to_graph as mol_to_graph_herg
from models.featurisation_ppb import mol_to_graph as mol_to_graph_ppb
from models.model_defs import GINRegressor, GATv2Regressor, GCN, PPBFuModel

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def get_available_adme_models():
    model_dir = Path(__file__).resolve().parent.parent / "models"
    return {
        "cyp3a4": {
            "path": model_dir / "cyp_3a4_gin.pt",
            "class": GINRegressor,
            "featuriser": mol_to_graph_cyp
        },
        "caco2": {
            "path": model_dir / "caco2_abpapp_gin.pt",
            "class": GINRegressor,
            "featuriser": mol_to_graph_cyp
        },
        "logd": {
            "path": model_dir / "logd_gat.pt",
            "class": GATv2Regressor,
            "featuriser": mol_to_graph_logd
        },
        "herg": {
            "path": model_dir / "herg_gnn_latest.pt",
            "class": GCN,
            "apply_sigmoid": True,
            "featuriser": mol_to_graph_herg
        },
            "ppb_f_u": {
            "path": model_dir / "ppb_f_u_gin.pt",
            "class": PPBFuModel,
            "featuriser": mol_to_graph_ppb
        },
    }

@register_task("adme_prediction", category="Prediction", description="Predict ADME properties (hERG, LogD, CYP3A4, A->B (Papp)) from SMILES using PyTorch models")
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

    if data is not None and "df" in data:
        df = data["df"].copy()
    else:
        df = pd.read_csv(input_file)

    if "smiles" not in df.columns:
        raise ValueError("Input file must contain a 'smiles' column.")

    df["smiles"] = df["smiles"].astype(str).str.strip()
    df = df[df["smiles"].notna() & df["smiles"].str.len() > 0].copy()

    result_df = df.copy()   # keep all columns

    device = torch.device("cpu" if torch.cuda.is_available() else "cpu")
    logger.debug(f"Using device: {device}")

    for model_name, info in model_configs.items():
        model_path = info["path"]
        model_class = info["class"]
        logger.debug(f"Loading model '{model_name}' from {model_path}")

        try:
            checkpoint = torch.load(model_path, map_location=device)

            required_keys = ['input_dim', 'model_state_dict']
            missing_keys = [k for k in required_keys if k not in checkpoint]
            if missing_keys:
                logger.error(f"Checkpoint missing keys: {missing_keys}")
                continue

            model_args = {
                "input_dim": checkpoint["input_dim"],
                "hidden_dim": checkpoint.get("hidden_dim", 128),
            }

            if model_class in [GINRegressor, GATv2Regressor]:
                model_args["global_feat_dim"] = checkpoint.get("global_feat_dim", 0)

            if model_class in [GATv2Regressor, GCN]:
                if "edge_dim" not in checkpoint:
                    logger.error(f"Missing 'edge_dim' in checkpoint for model '{model_name}'")
                    continue
                model_args["edge_dim"] = checkpoint["edge_dim"]

            model = model_class(**model_args)
            load_result = model.load_state_dict(checkpoint["model_state_dict"])
            logger.debug(f"load_state_dict result for '{model_name}': {load_result}")
            model.to(device)
            model.eval()

        except Exception as e:
            logger.exception(f"Error loading model '{model_name}': {e}")

        preds = []
        featuriser = info.get("featuriser") # get the featuriser for the specific model
        apply_sigmoid = info.get("apply_sigmoid", False) # flag for whether to apply sigmoind transform for 0/1 models

        for smi in df["smiles"]:
            try:
                graph_data = featuriser(smi)
                if graph_data is None:
                    logger.error(f"[{model_name}] Featuriser returned None for SMILES: '{smi}'")
                    preds.append(None)
                    continue
                
                for attr in ["x", "edge_index", "batch"]:
                    if not hasattr(data, attr):
                        logger.error(f"[{model_name}] Data object missing '{attr}' attribute for SMILES: {smi}")
                    elif getattr(data, attr) is None:
                        logger.error(f"[{model_name}] Data.{attr} is None for SMILES: {smi}")

                if hasattr(data, "x") and torch.isnan(data.x).any():
                    logger.warning(f"[{model_name}] data.x contains NaN for SMILES: {smi}")

                loader = DataLoader([graph_data], batch_size=1, shuffle=False)

                with torch.no_grad():
                    for batch in loader:
                        batch = batch.to(device)
                        output = model(batch)
                        if apply_sigmoid:
                            output = torch.sigmoid(output)
                        if output is None:
                            logger.error(f"[{model_name}] Model returned None for SMILES: {smi}")
                            preds.append(None)
                        else:
                            preds.append(output.item() if isinstance(output, torch.Tensor) else float(output))

            except Exception as e:
                logger.error(f"[{model_name}] Inference failed for SMILES '{smi}': {e}")
                preds.append(None)


        result_df[model_name] = preds

    return (Path(input_file).stem, result_df)
