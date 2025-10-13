import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GINConv, global_add_pool

class PPBFuModel(nn.Module):
    def __init__(self, input_dim, hidden_dim=64, global_feat_dim=0, n_layers=3, dropout=0.2):
        super().__init__()
        self.layers = nn.ModuleList()
        for i in range(n_layers):
            in_dim = input_dim if i == 0 else hidden_dim
            mlp = nn.Sequential(
                nn.Linear(in_dim, hidden_dim),
                nn.ReLU(),
                nn.Linear(hidden_dim, hidden_dim)
            )
            self.layers.append(GINConv(mlp))
        self.bns = nn.ModuleList([nn.BatchNorm1d(hidden_dim) for _ in range(n_layers)])

        self.readout = nn.Sequential(
            nn.Linear(hidden_dim + global_feat_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, 1)
        )

    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch
        for conv, bn in zip(self.layers, self.bns):
            x = F.relu(bn(conv(x, edge_index)))
        x = global_add_pool(x, batch)
        if hasattr(data, "global_features"):
            x = torch.cat([x, data.global_features.to(x.device)], dim=1)
        return self.readout(x).squeeze(1)
