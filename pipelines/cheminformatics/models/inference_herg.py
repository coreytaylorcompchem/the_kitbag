import torch
from torch_geometric.loader import DataLoader
from featurisation_herg import mol_to_graph
from model_defs import GCN
import pandas as pd

def load_model(model_path, input_dim, edge_dim, hidden_dim, device):
    checkpoint = torch.load('herg_gnn.pt', map_location=device)
    input_dim = checkpoint['input_dim']
    edge_dim = checkpoint['edge_dim']
    hidden_dim = checkpoint['hidden_dim']

    model = GCN(input_dim=input_dim, edge_dim=edge_dim, hidden_dim=hidden_dim)
    state = torch.load("herg_gnn.pt", map_location=device)
    model.load_state_dict(state['model_state_dict'])
    model.to(device)
    model.eval()
    return model

def predict_smiles_list(smiles_list, model, device, batch_size=32):
    data_list = [mol_to_graph(smi) for smi in smiles_list]
    data_list = [d for d in data_list if d is not None]

    loader = DataLoader(data_list, batch_size=batch_size, shuffle=False)
    preds = []

    with torch.no_grad():
        for batch in loader:
            batch = batch.to(device)
            out = model(batch)
            probs = torch.sigmoid(out)
            preds.extend(probs.cpu().numpy().tolist())

    return preds

def main(config):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    model = load_model(
        config["model_path"],
        config["input_dim"],
        config["edge_dim"],
        config["hidden_dim"],
        device
    )

    preds = predict_smiles_list(config["smiles_list"], model, device)
    df = pd.DataFrame({
        "smiles": config["smiles_list"],
        "prob_herg_block": preds
    })
    return df

if __name__ == "__main__":
    config = {
        "model_path": "herg_gnn_latest.pt",
        "input_dim": 30,
        "edge_dim": 6,
        "hidden_dim": 128,
        "smiles_list": ["CCO", "c1ccccc1", "CC(=O)O"]
    }

    df = main(config)
