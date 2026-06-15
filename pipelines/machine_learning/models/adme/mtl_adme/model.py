from torch_geometric.nn import GINEConv, global_add_pool, global_mean_pool
from torch_geometric.nn import GlobalAttention
from torch_geometric.nn import AttentionalAggregation
import torch.nn as nn
import torch.nn.functional as F
import torch

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

class GINRegressor(nn.Module):
    def __init__(self, input_dim, edge_dim, global_feat_dim, fp_dim, num_tasks,
                 hidden_dim=256, fp_hidden_dim=256, num_layers=5, dropout=0.1, task_groups=None, 
                 group_architecture=None, loss_name="huber",huber_delta=1.0,):
        super().__init__()
        self._debug_printed = False
        self.num_layers = num_layers
        self.dropout = dropout
        self.num_tasks = num_tasks 
        self.input_proj = nn.Linear(input_dim, hidden_dim)
        self.loss_name = loss_name
        self.huber_delta = huber_delta
        self.att_pool = AttentionalAggregation(
            gate_nn=nn.Sequential(
                nn.Linear(hidden_dim, hidden_dim // 2),
                nn.ReLU(),
                nn.Linear(hidden_dim // 2, 1)
            )
        )
        self.film = nn.Sequential(
            nn.Linear(global_feat_dim, hidden_dim * 2)
        )
        self.graph_norm = nn.LayerNorm(hidden_dim)
        self.fp_norm = nn.LayerNorm(fp_hidden_dim)
        self.global_norm = nn.LayerNorm(hidden_dim)
        self.group_architecture = (
            group_architecture or {}
        ) # store stage-specific parameters

        # =========================
        # TASK GROUPING
        # =========================
        if task_groups is None:
            raise ValueError("task_groups must be provided")

        self.task_groups = task_groups
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
            self.norms.append(nn.LayerNorm(hidden_dim))

        # =========================
        # Fingerprint MLP
        # =========================
        # self.fp_mlp = nn.Sequential(
        #     nn.Linear(fp_dim, 1024),
        #     nn.ReLU(),
        #     nn.Dropout(dropout), 
        #     nn.Linear(1024, fp_hidden_dim),
        #     nn.ReLU()
        # )
        self.fp_mlp = nn.Sequential(
            nn.Linear(fp_dim, 512),
            nn.LayerNorm(512),
            nn.ReLU(),
            nn.Dropout(dropout),

            nn.Linear(512, fp_hidden_dim),
            nn.LayerNorm(fp_hidden_dim),
            nn.ReLU()
        )

        # self.fp_gate = nn.Parameter(torch.zeros(fp_hidden_dim))
        # self.global_gate = nn.Parameter(torch.ones(hidden_dim)) # delete if things crash

        self.fp_gate = nn.Parameter(torch.ones(fp_hidden_dim))
        self.global_gate = nn.Parameter(torch.ones(hidden_dim))

        # =========================
        # Global features MLP
        # =========================
        self.global_proj = nn.Sequential(
            nn.Linear(global_feat_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout), 
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU()
        )

        # self.global_gate = nn.Parameter(torch.zeros(hidden_dim))

        # =========================
        # SHARED TRUNK
        # =========================

        final_dim = hidden_dim + fp_hidden_dim + hidden_dim

        self.pre_shared_norm = nn.LayerNorm(final_dim)

        self.shared = nn.Sequential(
            nn.Linear(final_dim, 512),
            nn.LayerNorm(512),
            nn.ReLU(),
            nn.Dropout(dropout)
        )

        # =========================
        # GROUP-SPECIFIC TRUNKS
        # =========================
        # self.group_trunks = nn.ModuleDict({
        #     "physchem": nn.Sequential(
        #         nn.Linear(512, 256),
        #         nn.ReLU(),
        #         nn.Dropout(dropout)
        #     ),
        #     "adme": nn.Sequential(
        #         nn.Linear(512, 256),
        #         nn.ReLU(),
        #         nn.Dropout(dropout)
        #     ),
        #     "cyp": nn.Sequential(
        #         nn.Linear(512, 256),
        #         nn.ReLU(),
        #         nn.Dropout(dropout)
        #     ),
        #     "tox": nn.Sequential(
        #         nn.Linear(512, 256),
        #         nn.ReLU(),
        #         nn.Dropout(dropout)
        #     ),
        # })

        self.group_trunks = nn.ModuleDict()

        task_output_dims = {}

        for group_name, task_indices in self.task_groups.items():

            cfg = self.group_architecture.get(
                group_name,
                {}
            )

            hidden_dim_group = cfg.get(
                "hidden_dim",
                256
            )

            output_dim_group = cfg.get(
                "output_dim",
                256
            )

            self.group_trunks[group_name] = nn.Sequential(

                nn.Linear(
                    512,
                    hidden_dim_group
                ),
                nn.ReLU(),
                nn.Dropout(dropout),

                nn.Linear(
                    hidden_dim_group,
                    output_dim_group
                ),
                nn.ReLU(),
                nn.Dropout(dropout),
            )

            for task_idx in task_indices:
                task_output_dims[task_idx] = output_dim_group
            
        self.heads = nn.ModuleList()

        for task_idx in range(num_tasks):

            head_input_dim = task_output_dims[
                task_idx
            ]

            self.heads.append(

                nn.Sequential(

                    nn.Linear(
                        head_input_dim,
                        64
                    ),

                    nn.ReLU(),
                    nn.Dropout(dropout),

                    nn.Linear(
                        64,
                        1
                    )
                )
            )

        logger.debug("Group architectures:")

        for group_name, trunk in self.group_trunks.items():

            first = trunk[0]
            second = trunk[3]

            logger.info(
                f"{group_name}: "
                f"{first.in_features}->{first.out_features} "
                f"-> "
                f"{second.out_features}"
            )

        # for group_name in self.task_groups.keys():

        #     cfg = self.group_architecture.get(
        #         group_name,
        #         {}
        #     )

        #     hidden_dim_group = cfg.get(
        #         "hidden_dim",
        #         256
        #     )

        #     output_dim_group = cfg.get(
        #         "output_dim",
        #         256
        #     )

        # self.group_trunks[group_name] = nn.Sequential(

        #     nn.Linear(512, hidden_dim_group),
        #     nn.ReLU(),
        #     nn.Dropout(dropout),

        #     nn.Linear(
        #         hidden_dim_group,
        #         output_dim_group
        #     ),
        #     nn.ReLU(),
        #     nn.Dropout(dropout),
        # )

        # =========================
        # Task-specific heads
        # =========================
        # self.heads = nn.ModuleList([
        #     nn.Sequential(
        #         nn.Linear(256, 128),
        #         nn.ReLU(),
        #         nn.Dropout(dropout),
        #         nn.Linear(128, 1)
        #     )
        #     for _ in range(num_tasks)
        # ])

        # try a smaller model
        # self.heads = nn.ModuleList([
        #     nn.Sequential(
        #         nn.Linear(256, 64),
        #         nn.ReLU(),
        #         nn.Dropout(dropout),
        #         nn.Linear(64, 1)
        #     )
        #     for _ in range(num_tasks)
        # ])

        # =========================
        # Uncertainty weights
        # =========================
        # self.log_vars = nn.Parameter(torch.zeros(num_tasks))

    def forward(self, data):
        x, edge_index, edge_attr, batch = (
            data.x, data.edge_index, data.edge_attr, data.batch
        )

        # =========================
        # GNN forward + residuals
        # =========================
        x = self.input_proj(x)
        global_feat = data.global_features

        if global_feat.dim() == 3:
            global_feat = global_feat.squeeze(1)

        film_params = self.film(global_feat)

        gamma, beta = torch.chunk(film_params, 2, dim=-1)

        gamma = 1.0 + 0.05 * torch.tanh(gamma)
        beta = 0.05 * beta

        gamma = gamma[batch]
        beta = beta[batch]

        for conv, norm in zip(self.convs, self.norms):

            residual = x

            h = conv(x, edge_index, edge_attr)
            h = norm(h)

            # FiLM conditioning
            h = gamma * h + beta

            # h = F.relu(h)
            # h = F.dropout(h, p=self.dropout, training=self.training)

            # x = residual + h

            h = F.relu(h)
            x = residual + h
            x = F.dropout(x, p=self.dropout, training=self.training)

        # Pooling
        # x = global_mean_pool(x, batch)
        x = self.att_pool(x, batch)

        # =========================
        # Fingerprints
        # =========================
        fp = data.fp
        if fp.dim() == 3:
            fp = fp.squeeze(1)

        fp = fp / (fp.norm(p=2, dim=1, keepdim=True) + 1e-6)
        # fp_out = self.fp_mlp(fp)
        if self.training:
            fp_mask = (torch.rand_like(fp) > 0.15).float()
            fp = fp * fp_mask

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
        # x = torch.cat([
        #     x,
        #     fp_out * torch.sigmoid(self.fp_gate),
        #     global_out * torch.sigmoid(self.global_gate)
        # ], dim=1)

        x_graph = self.graph_norm(x)

        x_fp = self.fp_norm(
            fp_out * torch.sigmoid(self.fp_gate)
        )

        x_global = self.global_norm(
            global_out * torch.sigmoid(self.global_gate)
        )

        x = torch.cat([
            x_graph,
            x_fp,
            x_global
        ], dim=1)

        # x = F.layer_norm(x, x.shape[1:])
        x = self.pre_shared_norm(x)

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
            # logger.debug(f"[ROUTING] Group '{group_name}' → tasks {task_indices}")
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
                    logger.debug(f"[HEAD] Task {task_idx} uses group '{group_name}'")

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