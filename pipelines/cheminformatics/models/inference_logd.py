import torch
from torch_geometric.loader import DataLoader
import pandas as pd
from pathlib import Path

from models.model_defs import GATv2Regressor
from models.featurisation_logd import mol_to_graph, descriptor_functions
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=True, simple_format=True)

def load_model(model_path, device):
    logger.debug(f"Loading model from: {model_path} to device: {device}")
    checkpoint = torch.load(model_path, map_location=device)

    try:
        model = GATv2Regressor(
            input_dim=checkpoint["input_dim"],
            edge_dim=checkpoint["edge_dim"],
            hidden_dim=checkpoint["hidden_dim"],
            global_feat_dim=checkpoint["global_feat_dim"]
        )
    except Exception as e:
        logger.error(f"Failed to instantiate GATv2Regressor: {e}", exc_info=True)
        raise  # Optional: re-raise if you want it to stop here
    
    load_result = model.load_state_dict(checkpoint["model_state_dict"])
    logger.debug(f"load_state_dict result: {load_result}")
    
    model.to(device)
    model.eval()
    logger.debug(f"Model loaded and moved to {device}")

    return model

def predict_smiles_list(smiles_list, model, device, batch_size=32):
    """
    Converts SMILES list to graphs and performs batch inference.
    """
    data_list = []
    for smi in smiles_list:
        data = mol_to_graph(smi)
        if data is not None:
            data_list.append(data)

    if not data_list:
        raise ValueError("No valid molecules in the input SMILES list.")

    loader = DataLoader(data_list, batch_size=batch_size, shuffle=False)

    preds = []
    with torch.no_grad():
        for batch in loader:
            batch = batch.to(device)
            out = model(batch)
            preds.extend(out.cpu().numpy().tolist())

    return preds

def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # Path to the trained model checkpoint
    model_path = 'models/logd_gat.pt'

    # List of SMILES for prediction
    smiles_list = ["CCO", "c1ccccc1", "CC(=O)O"]

    # Load model and predict
    model = load_model(model_path, device)
    preds = predict_smiles_list(smiles_list, model, device)

    # Output results
    df = pd.DataFrame({'smiles': smiles_list, 'predicted_LogD': preds})
    print(df)
    return df

if __name__ == "__main__":
    main()
