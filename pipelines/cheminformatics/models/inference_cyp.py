import torch
from model_defs import GINRegressor
from featurization import mol_to_graph, descriptor_functions
from torch_geometric.loader import DataLoader
import pandas as pd
from pathlib import Path

def load_model(model_path, input_dim, global_feat_dim, device):
    model = GINRegressor(input_dim=input_dim, global_feat_dim=global_feat_dim)
    model.load_state_dict(torch.load(model_path, map_location=device))
    model.to(device)
    model.eval()
    return model

def predict_smiles_list(smiles_list, model, device, batch_size=32):
    data_list = []
    for smi in smiles_list:
        data = mol_to_graph(smi)
        if data is not None:
            data_list.append(data)

    loader = DataLoader(data_list, batch_size=batch_size, shuffle=False)

    preds = []
    with torch.no_grad():
        for batch in loader:
            batch = batch.to(device)
            out = model(batch)
            preds.extend(out.cpu().numpy().tolist())

    return preds

def main(config):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # Assumes you know the input_dim and global_feat_dim (can be read from a sample data)
    input_dim = config['input_dim']
    global_feat_dim = config['global_feat_dim']

    model = load_model(config['model_path'], input_dim, global_feat_dim, device)

    smiles_list = config['smiles_list']

    preds = predict_smiles_list(smiles_list, model, device)

    # Return or save predictions as needed, e.g. as dataframe
    df = pd.DataFrame({'smiles': smiles_list, 'predicted_value': preds})
    print(df)
    return df

if __name__ == "__main__":
    # Example usage:
    example_smiles = ["CCO", "c1ccccc1", "CC(=O)O"]
    config = {
        'model_path': 'models/best_model.pth',
        'input_dim': 60,  # example input_dim (adjust according to your featurization)
        'global_feat_dim': len(descriptor_functions),
        'smiles_list': example_smiles,
    }
    main(config)
