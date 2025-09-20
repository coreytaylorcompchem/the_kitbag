import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GINConv, GATv2Conv, global_mean_pool

class GINRegressor(nn.Module):
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

        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = self.conv2(x, edge_index)
        x = F.relu(x)

        x = global_mean_pool(x, batch)

        if hasattr(data, 'global_features'):
            global_feat = data.global_features.to(x.device)
            x = torch.cat([x, global_feat], dim=1)

        return self.lin(x).squeeze(1)

class GATv2Regressor(torch.nn.Module):
    def __init__(self, input_dim, edge_dim, hidden_dim=128, heads=2, global_feat_dim=0):
        # print(f"[DEBUG] Entering GATv2Regressor.__init__()")
        super().__init__()
        # print(f"[DEBUG] super().__init__() complete")
        self.edge_encoder = nn.Linear(edge_dim, hidden_dim)

        self.gat1 = GATv2Conv(input_dim, hidden_dim, heads=heads, dropout=0.1, edge_dim=hidden_dim)
        self.gat2 = GATv2Conv(hidden_dim * heads, hidden_dim, heads=1, concat=True, dropout=0.1, edge_dim=hidden_dim)

        self.global_mlp = nn.Sequential(
            nn.Linear(global_feat_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU()
        )

        self.lin = nn.Sequential(
            nn.Linear(hidden_dim + (hidden_dim // 2), hidden_dim // 2),
            nn.ReLU(),
            nn.Linear(hidden_dim // 2, 1)
        )

        self.attentions = {}  # Initialize attention storage

    def forward(self, data):
        x, edge_index, edge_attr, batch = data.x, data.edge_index, data.edge_attr, data.batch
        device = next(self.parameters()).device

        x, edge_index, edge_attr, batch = x.to(device), edge_index.to(device), edge_attr.to(device), batch.to(device)
        edge_attr = self.edge_encoder(edge_attr)

        # Extract attention weights from GATv2Conv layers
        x, attn1 = self.gat1(x, edge_index, edge_attr, return_attention_weights=True)
        self.attentions['layer1'] = attn1  # attn1 is tuple (edge_index, attention_weights)
        x = F.elu(x)

        x, attn2 = self.gat2(x, edge_index, edge_attr, return_attention_weights=True)
        self.attentions['layer2'] = attn2
        x = F.elu(x)

        x = global_mean_pool(x, batch)  # [batch_size, hidden_dim]

        if hasattr(data, 'global_features'):
            global_features = data.global_features.to(x.device)
            g_feat = self.global_mlp(global_features)
            x = torch.cat([x, g_feat], dim=1)

        return self.lin(x).squeeze(1)

    
    def forward_with_attention(self, data):
        x, edge_index, edge_attr, batch = data.x, data.edge_index, data.edge_attr, data.batch
        device = next(self.parameters()).device

        x, edge_index, edge_attr, batch = x.to(device), edge_index.to(device), edge_attr.to(device), batch.to(device)
        edge_attr = self.edge_encoder(edge_attr)

        # Capture attention weights here (GATv2Conv returns them if return_attention_weights=True)
        # Layer 1 with attention extraction
        x, (edge_idx1, attn1) = self.gat1(x, edge_index, edge_attr, return_attention_weights=True)
        self.attentions['layer1'] = attn1  # store attention (tuple: (edge_index, attention_weights))
        x = F.elu(x)

        # Layer 2 with attention extraction
        x, attn2 = self.gat2(x, edge_index, edge_attr, return_attention_weights=True)
        self.attentions['layer2'] = attn2
        x = F.elu(x)

        x = global_mean_pool(x, batch)

        if hasattr(data, 'global_features'):
            global_features = data.global_features.to(x.device)
            g_feat = self.global_mlp(global_features)
            x = torch.cat([x, g_feat], dim=1)

        out = self.lin(x).squeeze(1)
        return out, attn1, edge_idx1