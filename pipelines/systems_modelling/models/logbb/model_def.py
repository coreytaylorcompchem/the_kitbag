import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GINEConv, global_add_pool


class LogBBModel(nn.Module):
    """GIN model for logBB (brain–blood partition coefficient) prediction."""

    def __init__(self, input_dim, edge_dim, global_feat_dim,
                 hidden_dim=128, n_layers=4, dropout=0.2):
        super().__init__()

        self.layers = nn.ModuleList()
        self.bns = nn.ModuleList()

        for i in range(n_layers):
            in_dim = input_dim if i == 0 else hidden_dim
            mlp = nn.Sequential(
                nn.Linear(in_dim, hidden_dim),
                nn.LeakyReLU(0.1),
                nn.Linear(hidden_dim, hidden_dim)
            )
            self.layers.append(GINEConv(mlp, edge_dim=edge_dim))
            self.bns.append(nn.BatchNorm1d(hidden_dim))

        # Always define projection from input_dim → hidden_dim
        # used when shapes differ (runtime check)
        self.res_proj = nn.Linear(input_dim, hidden_dim) if input_dim != hidden_dim else nn.Identity()

        # Global descriptor projection
        self.global_proj = nn.Sequential(
            nn.Linear(global_feat_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout)
        )

        # Readout head
        self.readout = nn.Sequential(
            nn.Linear(hidden_dim * 2, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Linear(hidden_dim // 2, 1)
        )

    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch
        edge_attr = getattr(data, "edge_attr", None)

        for conv, bn in zip(self.layers, self.bns):
            x_res = x
            x = F.leaky_relu(bn(conv(x, edge_index, edge_attr)))

            # Always project when feature dims mismatch
            if x_res.shape[1] != x.shape[1]:
                x_res = self.res_proj(x_res)

            x = x + x_res  # residual

        # Global pooling
        x = global_add_pool(x, batch)
        x = F.layer_norm(x, x.shape[-1:])

        # Global molecular features
        if hasattr(data, "global_features"):
            gf = self.global_proj(data.global_features.to(x.device))
            x = torch.cat([x, gf], dim=1)

        return self.readout(x).squeeze(1)

