from torch_geometric.nn import GINConv, global_mean_pool
import torch.nn as nn
import torch.nn.functional as F
import torch

from torch_geometric.nn import GINEConv, global_add_pool, BatchNorm

class GINRegressor(nn.Module):
    def __init__(self, input_dim, edge_dim, global_feat_dim, fp_dim, num_tasks,
                 hidden_dim=512, fp_hidden_dim=512, num_layers=5, dropout=0.05):
        super().__init__()
        self.num_layers = num_layers
        self.dropout = dropout

        # GIN layers
        self.convs = nn.ModuleList()
        self.norms = nn.ModuleList()
        for i in range(num_layers):
            in_dim = input_dim if i==0 else hidden_dim
            nn_layer = nn.Sequential(
                nn.Linear(in_dim, hidden_dim),
                nn.ReLU(),
                nn.Linear(hidden_dim, hidden_dim)
            )
            self.convs.append(GINEConv(nn_layer, edge_dim=edge_dim))
            self.norms.append(BatchNorm(hidden_dim))

        # Fingerprint MLP 
        self.fp_mlp = nn.Sequential(
            nn.Linear(fp_dim, 1024),
            nn.ReLU(),
            nn.Linear(1024, 512),
            nn.ReLU(),
            nn.Linear(512, fp_hidden_dim),
            nn.ReLU()
        )

        # Global features MLP 
        self.global_proj = nn.Sequential(
            nn.Linear(global_feat_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU()
        )

        # Final predictor 
        final_dim = hidden_dim + fp_hidden_dim + hidden_dim
        self.lin = nn.Sequential(
            nn.Linear(final_dim, 512),
            nn.ReLU(),
            nn.Linear(512, num_tasks)
        )

    def forward(self, data):
        x, edge_index, edge_attr, batch = data.x, data.edge_index, data.edge_attr, data.batch

        # GIN layers 
        for conv, norm in zip(self.convs, self.norms):
            x = conv(x, edge_index, edge_attr)
            x = norm(x)
            x = F.relu(x)

        # Sum pooling
        x = global_add_pool(x, batch)

        # Fingerprint 
        fp = data.fp
        if fp.dim() == 3:
            fp = fp.squeeze(1)
        fp = fp / (fp.norm(p=2, dim=1, keepdim=True)+1e-6)
        fp_out = self.fp_mlp(fp)

        # Global features 
        global_feat = data.global_features
        if global_feat.dim() == 3:
            global_feat = global_feat.squeeze(1)
        global_out = self.global_proj(global_feat)

        # Concatenate 
        x = torch.cat([x, fp_out, global_out], dim=1)
        return self.lin(x)   # shape (batch, num_tasks)