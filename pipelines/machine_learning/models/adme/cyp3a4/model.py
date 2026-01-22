from torch_geometric.nn import GINConv, global_mean_pool
import torch.nn as nn
import torch.nn.functional as F
import torch

class CYP3A4GIN(nn.Module):
    def __init__(self, input_dim, hidden_dim, global_feat_dim):
        super().__init__()

        nn1 = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim)
        )
        nn2 = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim)
        )

        self.conv1 = GINConv(nn1)
        self.conv2 = GINConv(nn2)

        self.readout = nn.Sequential(
            nn.Linear(hidden_dim + global_feat_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 1)
        )

    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch
        x = F.relu(self.conv1(x, edge_index))
        x = F.relu(self.conv2(x, edge_index))
        x = global_mean_pool(x, batch)

        if hasattr(data, "global_features"):
            x = torch.cat([x, data.global_features.to(x.device)], dim=1)

        return self.readout(x).squeeze(1)
