from torch_geometric.nn import GINEConv, global_add_pool, BatchNorm, global_mean_pool
import torch.nn as nn
import torch.nn.functional as F
import torch

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

class GINRegressor(nn.Module):
    def __init__(self, input_dim, edge_dim, global_feat_dim, fp_dim, num_tasks,
                 hidden_dim=512, fp_hidden_dim=512, num_layers=5, dropout=0.1):
        super().__init__()
        self._debug_printed = False
        self.num_layers = num_layers
        self.dropout = dropout
        self.num_tasks = num_tasks 
        self.input_proj = nn.Linear(input_dim, hidden_dim)

        # =========================
        # TASK GROUPING (for grouping architecture)
        # =========================
        self.task_groups = {
            "physchem": [0, 1, 2],
            "adme": [3, 4],
            "cyp": [5, 6, 7, 8],
            "tox": [9],
        }
        # =========================
        # DEBUG: validate task grouping
        # =========================
        all_grouped_tasks = sorted([
            t for group in self.task_groups.values() for t in group
        ])

        if all_grouped_tasks != list(range(self.num_tasks)):
            raise ValueError(
                f"Task grouping mismatch.\n"
                f"Expected: {list(range(self.num_tasks))}\n"
                f"Got: {all_grouped_tasks}"
            )

        logger.debug("Task group configuration:")
        for group, tasks in self.task_groups.items():
            logger.debug(f"  {group}: tasks {tasks}")

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
        # final_dim = hidden_dim + fp_hidden_dim + hidden_dim

        # self.shared = nn.Sequential(
        #     nn.Linear(final_dim, 512),
        #     nn.ReLU(),
        #     nn.Dropout(dropout)
        # )

        # =========================
        # SHARED TRUNK
        # =========================

        final_dim = hidden_dim + fp_hidden_dim + hidden_dim

        self.shared = nn.Sequential(
            nn.Linear(final_dim, 512),
            nn.ReLU(),
            nn.Dropout(dropout)
        )

        # =========================
        # GROUP-SPECIFIC TRUNKS
        # =========================
        self.group_trunks = nn.ModuleDict({
            "physchem": nn.Sequential(
                nn.Linear(512, 256),
                nn.ReLU(),
                nn.Dropout(dropout)
            ),
            "adme": nn.Sequential(
                nn.Linear(512, 256),
                nn.ReLU(),
                nn.Dropout(dropout)
            ),
            "cyp": nn.Sequential(
                nn.Linear(512, 256),
                nn.ReLU(),
                nn.Dropout(dropout)
            ),
            "tox": nn.Sequential(
                nn.Linear(512, 256),
                nn.ReLU(),
                nn.Dropout(dropout)
            ),
        })

        # =========================
        # Task-specific heads
        # =========================
        self.heads = nn.ModuleList([
            nn.Sequential(
                # nn.Linear(512, 128),
                nn.Linear(256, 128),
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
        x = self.input_proj(x)

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
        x = torch.cat([
            x,
            fp_out * torch.sigmoid(self.fp_gate),
            global_out * torch.sigmoid(self.global_gate)
        ], dim=1)

        # =========================
        # Shared trunk
        # =========================
        x_shared = self.shared(x)

        # =========================
        # DEBUG BLOCK 1: shared representation
        # =========================
        if not self._debug_printed:
            logger.debug("Shared representation:")
            logger.debug(f"  shape: {x_shared.shape}")
            logger.debug(f"  mean: {x_shared.mean().item():.4f}")
            logger.debug(f"  std:  {x_shared.std().item():.4f}")

        # Prepare output container
        outputs = [None] * self.num_tasks

        # =========================
        # Group routing
        # =========================
        for group_name, task_indices in self.task_groups.items():
            x_group = self.group_trunks[group_name](x_shared)

            # =========================
            # DEBUG BLOCK 2: per-group representation
            # =========================
            if not self._debug_printed:
                logger.debug(f"Group '{group_name}' representation:")
                logger.debug(f"  tasks: {task_indices}")
                logger.debug(f"  shape: {x_group.shape}")
                logger.debug(f"  mean: {x_group.mean().item():.4f}")
                logger.debug(f"  std:  {x_group.std().item():.4f}")

            for task_idx in task_indices:
                out_i = self.heads[task_idx](x_group)
                outputs[task_idx] = out_i

                # =========================
                # DEBUG BLOCK 3: per-task output
                # =========================
                if not self._debug_printed:
                    logger.debug(f"Task {task_idx} output mean: {out_i.mean().item():.4f}")

        # Final output
        out = torch.cat(outputs, dim=1)

        # =========================
        # DEBUG BLOCK 4: final output
        # =========================
        if not self._debug_printed:
            logger.debug(f"Final output shape: {out.shape}")
            logger.debug(f"Output mean: {out.mean().item():.4f}")
            logger.debug(f"Output std: {out.std().item():.4f}")

            self._debug_printed = True

        return out