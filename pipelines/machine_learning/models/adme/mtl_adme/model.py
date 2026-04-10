from torch_geometric.nn import GINEConv, global_add_pool, BatchNorm, global_mean_pool
import torch.nn as nn
import torch.nn.functional as F
import torch


class GINRegressor(nn.Module):
    def __init__(self, input_dim, edge_dim, global_feat_dim, fp_dim, num_tasks,
                 hidden_dim=512, fp_hidden_dim=512, num_layers=5, dropout=0.1):
        super().__init__()
        self.num_layers = num_layers
        self.dropout = dropout
        self.num_tasks = num_tasks 
        self.input_proj = nn.Linear(input_dim, hidden_dim)

        # =========================
        # GNN backbone
        # =========================
        self.convs = nn.ModuleList()
        self.norms = nn.ModuleList()

        for _ in range(num_layers):
            nn_layer = nn.Sequential(
                nn.Linear(hidden_dim, hidden_dim),
                nn.ReLU(),
                nn.Linear(hidden_dim, hidden_dim)
            )
            self.convs.append(GINEConv(nn_layer, edge_dim=edge_dim))
            self.norms.append(BatchNorm(hidden_dim))

        # =========================
        # Fingerprint MLP
        # =========================
        self.fp_mlp = nn.Sequential(
            nn.Linear(fp_dim, 1024),
            nn.ReLU(),
            nn.Dropout(dropout), 
            nn.Linear(1024, fp_hidden_dim),
            nn.ReLU()
        )

        self.fp_gate = nn.Parameter(torch.zeros(fp_hidden_dim))

        # =========================
        # Global features MLP
        # =========================
        self.global_proj = nn.Sequential(
            nn.Linear(global_feat_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),  # NEW
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU()
        )

        self.global_gate = nn.Parameter(torch.zeros(hidden_dim))

        # =========================
        # Shared trunk
        # =========================
        final_dim = hidden_dim + fp_hidden_dim + hidden_dim

        self.shared = nn.Sequential(
            nn.Linear(final_dim, 512),
            nn.ReLU(),
            nn.Dropout(dropout)
        )

        # =========================
        # Task-specific heads
        # =========================
        self.heads = nn.ModuleList([
            nn.Sequential(
                nn.Linear(512, 128),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(128, 1)
            )
            for _ in range(num_tasks)
        ])

        # =========================
        # Uncertainty weights
        # =========================
        self.log_vars = nn.Parameter(torch.zeros(num_tasks))

    def forward(self, data):
        x, edge_index, edge_attr, batch = (
            data.x, data.edge_index, data.edge_attr, data.batch
        )

        # =========================
        # GNN forward + residuals
        # =========================
        # Project input once
        x = self.input_proj(x)

        # GNN forward + residuals
        for conv, norm in zip(self.convs, self.norms):
            h = conv(x, edge_index, edge_attr)
            h = norm(h)
            h = F.relu(h)
            x = x + h

        # Pooling
        x = global_mean_pool(x, batch)

        # =========================
        # Fingerprints
        # =========================
        fp = data.fp
        if fp.dim() == 3:
            fp = fp.squeeze(1)

        fp = fp / (fp.norm(p=2, dim=1, keepdim=True) + 1e-6)
        fp_out = self.fp_mlp(fp)

        # =========================
        # Global features
        # =========================
        global_feat = data.global_features
        if global_feat.dim() == 3:
            global_feat = global_feat.squeeze(1)

        global_out = self.global_proj(global_feat)

        # =========================
        # Combine
        # =========================
        alpha = torch.sigmoid(self.fp_gate)
        x = torch.cat([
            x,
            fp_out * torch.sigmoid(self.fp_gate),
            global_out * torch.sigmoid(self.global_gate)
        ], dim=1)

        # =========================
        # Shared trunk
        # =========================
        x = self.shared(x)

        # =========================
        # Task heads
        # =========================
        outputs = []
        for head in self.heads:
            outputs.append(head(x))

        out = torch.cat(outputs, dim=1)

        return out  # shape (batch, num_tasks)