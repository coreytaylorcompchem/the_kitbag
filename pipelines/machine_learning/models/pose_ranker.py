import torch
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, global_mean_pool


class PoseScoringGNN(torch.nn.Module):
    """
    Graph neural network for molecular pose scoring.

    Args:
        in_dim (int): input feature dimension (node features)
        hidden_dim (int): hidden layer size (default: 128)
        n_outputs (int): number of regression targets (default: 9)
    """
    def __init__(self, in_dim, hidden_dim=128, n_outputs=9):
        super().__init__()
        self.conv1 = GCNConv(in_dim, hidden_dim)
        self.conv2 = GCNConv(hidden_dim, hidden_dim)
        self.conv3 = GCNConv(hidden_dim, hidden_dim)
        self.lin1 = torch.nn.Linear(hidden_dim, hidden_dim // 2)
        self.lin2 = torch.nn.Linear(hidden_dim // 2, n_outputs)

    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch

        x = F.relu(self.conv1(x, edge_index))
        x = F.relu(self.conv2(x, edge_index))
        x = F.relu(self.conv3(x, edge_index))

        x = global_mean_pool(x, batch)
        x = F.relu(self.lin1(x))
        out = self.lin2(x)
        return out
