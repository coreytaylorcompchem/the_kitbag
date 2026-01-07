import torch
from torch_geometric.loader import DataLoader
from pathlib import Path

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

from models.featurisation_cyp import mol_to_graph as mol_to_graph_cyp
from models.featurisation_logd import mol_to_graph as mol_to_graph_logd
from models.featurisation_herg import mol_to_graph as mol_to_graph_herg
from models.model_defs import GINRegressor, GATv2Regressor, GCN

logger = setup_logger(__name__, debug_mode=False, simple_format=True)


# -------------------------
# Model registry
# -------------------------
def get_available_adme_models():
    model_dir = Path(__file__).resolve().parent.parent / "models"
    return {
        "cyp3a4": {
            "path": model_dir / "cyp_3a4_gin.pt",
            "class": GINRegressor,
            "featuriser": mol_to_graph_cyp,
        },
        "caco2": {
            "path": model_dir / "caco2_abpapp_gin.pt",
            "class": GINRegressor,
            "featuriser": mol_to_graph_cyp,
        },
        "logd": {
            "path": model_dir / "logd_gat.pt",
            "class": GATv2Regressor,
            "featuriser": mol_to_graph_logd,
        },
        "herg": {
            "path": model_dir / "herg_gnn_latest.pt",
            "class": GCN,
            "apply_sigmoid": True,
            "featuriser": mol_to_graph_herg,
        },
    }


# -------------------------
# Cached model loader
# -------------------------
def get_adme_models(backend, requested_models):
    if "adme_models" in backend.cache:
        return backend.cache["adme_models"]

    all_models = get_available_adme_models()
    device = torch.device("cpu")

    if not requested_models:
        requested_models = list(all_models.keys())

    models = {}
    for name in requested_models:
        info = all_models[name]
        checkpoint = torch.load(info["path"], map_location=device)

        model_args = {
            "input_dim": checkpoint["input_dim"],
            "hidden_dim": checkpoint.get("hidden_dim", 128),
        }

        if info["class"] in [GINRegressor, GATv2Regressor]:
            model_args["global_feat_dim"] = checkpoint.get("global_feat_dim", 0)

        if info["class"] in [GATv2Regressor, GCN]:
            model_args["edge_dim"] = checkpoint["edge_dim"]

        model = info["class"](**model_args)
        model.load_state_dict(checkpoint["model_state_dict"])
        model.to(device)
        model.eval()

        models[name] = {
            "model": model,
            "featuriser": info["featuriser"],
            "apply_sigmoid": info.get("apply_sigmoid", False),
            "device": device,
        }

        logger.debug(f"[ADME] Loaded model: {name}")

    backend.cache["adme_models"] = models
    return models


# -------------------------
# Workflow task
# -------------------------
@register_task(
    "adme_prediction",
    category="Prediction",
    description="Predict ADME properties per ligand (cached GNN inference).",
)
def adme_prediction(backend, ligand, config, **kwargs):

    smiles = ligand.get("smiles")
    if not smiles:
        raise ValueError(f"Ligand '{ligand.get('name')}' missing SMILES.")

    requested_models = config.get("adme", {}).get("models", [])
    models = get_adme_models(backend, requested_models)

    results = {}
    for name, m in models.items():
        try:
            graph = m["featuriser"](smiles)
            if graph is None:
                results[name] = None
                continue

            loader = DataLoader([graph], batch_size=1)
            with torch.no_grad():
                for batch in loader:
                    batch = batch.to(m["device"])
                    out = m["model"](batch)
                    if m["apply_sigmoid"]:
                        out = torch.sigmoid(out)
                    results[name] = float(out.item())

        except Exception as e:
            logger.error(f"[ADME:{name}] Failed for {ligand['name']}: {e}")
            results[name] = None

    ligand.setdefault("adme", {}).update(results)
