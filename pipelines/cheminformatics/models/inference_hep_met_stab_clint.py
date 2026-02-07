import torch
import pandas as pd
from pathlib import Path
from torch_geometric.loader import DataLoader

from model_defs import ClintModel
from featurisation_hep_met_stab_clint import (
    mol_to_graph,
    descriptor_functions,
    compiled_smarts,
)

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

# ==========================================================
# Model loading
# ==========================================================
def load_model(
    model_path: str,
    input_dim: int,
    global_feat_dim: int,
    hidden_dim: int,
    device: torch.device,
):
    """
    Load trained hepatic intrinsic clearance (Clint) model and print parameter shapes.
    """
    model = ClintModel(
        input_dim=input_dim,
        hidden_dim=hidden_dim,
        global_feat_dim=global_feat_dim,
    )

    state_dict = torch.load(
        model_path,
        map_location=device,
        weights_only=True,
    )

    print("Checkpoint parameter shapes:")
    for k, v in state_dict.items():
        print(f"{k}: {tuple(v.shape)}")

    try:
        model.load_state_dict(state_dict)
    except RuntimeError as e:
        print("Error loading state_dict:", e)

    print("\nModel's current parameter shapes:")
    for name, param in model.named_parameters():
        print(f"{name}: {tuple(param.shape)}")

    model.to(device)
    model.eval()

    print(f"\nLoaded hep-Clint model from {model_path}")
    return model



# ==========================================================
# Prediction
# ==========================================================
def predict_smiles_list(
    smiles_list,
    model,
    device,
    batch_size: int = 32,
):
    """
    Run hepatic Clint inference on a list of SMILES.
    """
    data_list = []
    valid_smiles = []

    for smi in smiles_list:
        data = mol_to_graph(smi)
        if data is None:
            logger.warning(f"Skipping invalid SMILES: {smi}")
            continue
        data_list.append(data)
        valid_smiles.append(smi)

    if not data_list:
        return [], []

    loader = DataLoader(
        data_list,
        batch_size=batch_size,
        shuffle=False,
    )

    preds = []
    with torch.no_grad():
        for batch in loader:
            batch = batch.to(device)
            out = model(batch)
            preds.extend(out.cpu().numpy().tolist())

    return valid_smiles, preds


# ==========================================================
# Main entry point
# ==========================================================
def main(config: dict):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    logger.info(f"Using device: {device}")

    input_dim = config["input_dim"]
    global_feat_dim = config.get(
        "global_feat_dim",
        len(descriptor_functions) + len(compiled_smarts),
    )

    model = load_model(
        model_path=config["model_path"],
        input_dim=input_dim,
        global_feat_dim=global_feat_dim,
        hidden_dim=config["hidden_dim"],
        device=device,
    )

    smiles_list = config["smiles_list"]

    valid_smiles, preds = predict_smiles_list(
        smiles_list=smiles_list,
        model=model,
        device=device,
        batch_size=config.get("batch_size", 32),
    )

    df = pd.DataFrame(
        {
            "smiles": valid_smiles,
            "hep_clint_pred": preds,
        }
    )

    if "output_path" in config:
        output_path = Path(config["output_path"])
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path, index=False)
        logger.info(f"Saved predictions to {output_path}")

    return df


# ==========================================================
# Example usage
# ==========================================================
if __name__ == "__main__":
    example_smiles = [
        "CCO",
        "c1ccccc1",
        "CC(=O)O",
    ]

    config = {
        "model_path": "models/hep_clint_best_model.pt",
        "input_dim": 44,     # MUST match training atom feature size
        "hidden_dim": 128,    # or 128 — must match checkpoint
        "global_feat_dim": 16,
        "smiles_list": example_smiles,
        "batch_size": 32,
        # "output_path": "outputs/hep_clint_predictions.csv",
    }

    df = main(config)
    
    # print(df.head())
