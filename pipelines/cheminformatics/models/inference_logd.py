import torch
from torch_geometric.loader import DataLoader
import pandas as pd
from pathlib import Path

from model_defs import GATv2Regressor
from featurisation_logd import mol_to_graph, descriptor_functions

def load_model(model_path, device):
    """
    Load model from checkpoint dictionary.
    """
    checkpoint = torch.load(model_path, map_location=device)

    model = GATv2Regressor(
        input_dim=checkpoint["input_dim"],
        edge_dim=checkpoint["edge_dim"],
        hidden_dim=checkpoint["hidden_dim"],
        global_feat_dim=checkpoint["global_feat_dim"]
    )

    model.load_state_dict(checkpoint["model_state_dict"])
    model.to(device)
    model.eval()

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
