from torch_geometric.nn import GINConv, global_mean_pool
import torch.nn as nn
import torch.nn.functional as F
import torch

# Trying a lighter model than an attention-based model

class GINRegressor(torch.nn.Module):
    def __init__(self, input_dim, hidden_dim=64, global_feat_dim=0):
        super().__init__()

        nn1 = nn.Sequential(nn.Linear(input_dim, hidden_dim), nn.ReLU(), nn.Linear(hidden_dim, hidden_dim))
        nn2 = nn.Sequential(nn.Linear(hidden_dim, hidden_dim), nn.ReLU(), nn.Linear(hidden_dim, hidden_dim))

        self.conv1 = GINConv(nn1)
        self.conv2 = GINConv(nn2)

        self.lin = nn.Sequential(
            nn.Linear(hidden_dim + global_feat_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 1)
        )

    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch

        # Ensure x, edge_index, batch are on same device
        x = x.to(batch.device)
        edge_index = edge_index.to(batch.device)
        batch = batch.to(batch.device)

        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = self.conv2(x, edge_index)
        x = F.relu(x)

        x = global_mean_pool(x, batch)

        if hasattr(data, 'global_features'):
            # Repeat global features per node batch
            global_feat = data.global_features.to(x.device).float()
            if global_feat.dim() == 1:
                # shape [num_global_features] -> [batch_size, num_global_features]
                global_feat = global_feat.unsqueeze(0).repeat(x.size(0), 1)
            x = torch.cat([x, global_feat], dim=1)

        return self.lin(x).squeeze(1)
