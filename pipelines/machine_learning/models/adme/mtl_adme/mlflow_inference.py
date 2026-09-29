from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch

from torch_geometric.data import Batch
from torch_geometric.loader import DataLoader

from models.adme.mtl_adme.featurisation import (
    mol_to_graph,
    prepare_features,
)

from models.adme.mtl_adme.model import (
    GINRegressor,
)

from modules.utils.transforms import (
    log_to_vd,
    logfu_to_ppb,
    logit_to_bioavailability,
    pic50_to_ic50,
)

REQUIRED_CHECKPOINT_KEYS = {
    "model_state_dict",
    "input_dim",
    "edge_dim",
    "global_feat_dim",
    "fp_dim",
    "num_tasks",
    "params",
    "task_names",
    "task_groups",
    "label_scalers",
    "label_transform_metadata",
    "global_feature_scaler",
}


def load_checkpoint(
    checkpoint_path: str | Path,
) -> dict[str, Any]:
    """
    Load and validate a trusted ADME checkpoint.

    The checkpoint contains sklearn preprocessing objects, so
    weights_only=False is required. Only load checkpoints created by
    the trusted training pipeline.
    """
    checkpoint_path = Path(
        checkpoint_path
    ).expanduser().resolve()

    if not checkpoint_path.is_file():
        raise FileNotFoundError(
            f"Checkpoint does not exist: {checkpoint_path}"
        )

    checkpoint = torch.load(
        checkpoint_path,
        map_location="cpu",
        weights_only=False,
    )

    if not isinstance(checkpoint, dict):
        raise TypeError(
            "Expected the ADME checkpoint to contain a dictionary."
        )

    missing_keys = (
        REQUIRED_CHECKPOINT_KEYS
        - set(checkpoint.keys())
    )

    if missing_keys:
        raise KeyError(
            "Checkpoint is missing required inference metadata: "
            f"{sorted(missing_keys)}"
        )

    task_names = list(
        checkpoint["task_names"]
    )

    num_tasks = int(
        checkpoint["num_tasks"]
    )

    if len(task_names) != num_tasks:
        raise ValueError(
            f"Checkpoint contains {len(task_names)} "
            f"task names but num_tasks={num_tasks}."
        )

    label_scalers = checkpoint[
        "label_scalers"
    ]

    transform_metadata = checkpoint[
        "label_transform_metadata"
    ]

    missing_scalers = [
        task_name
        for task_name in task_names
        if task_name not in label_scalers
    ]

    if missing_scalers:
        raise KeyError(
            "Checkpoint is missing label scalers "
            f"for tasks: {missing_scalers}"
        )

    missing_transforms = [
        task_name
        for task_name in task_names
        if task_name not in transform_metadata
    ]

    if missing_transforms:
        raise KeyError(
            "Checkpoint is missing transform metadata "
            f"for tasks: {missing_transforms}"
        )

    if checkpoint["global_feature_scaler"] is None:
        raise ValueError(
            "Checkpoint contains "
            "global_feature_scaler=None. "
            "Inference requires the fitted scaler "
            "used during training."
        )  

    return checkpoint

def build_model_from_checkpoint(
    checkpoint: dict[str, Any],
    device: torch.device,
) -> GINRegressor:
    """
    Reconstruct GINRegressor exactly from checkpoint metadata.
    """
    params = dict(
        checkpoint["params"]
    )

    # Training-only parameters must not reach the model constructor.
    params.pop("lr", None)
    params.pop("weight_decay", None)

    model = GINRegressor(
        input_dim=int(
            checkpoint["input_dim"]
        ),
        edge_dim=int(
            checkpoint["edge_dim"]
        ),
        global_feat_dim=int(
            checkpoint["global_feat_dim"]
        ),
        fp_dim=int(
            checkpoint["fp_dim"]
        ),
        num_tasks=int(
            checkpoint["num_tasks"]
        ),
        task_groups=checkpoint[
            "task_groups"
        ],
        group_architecture=checkpoint.get(
            "group_architecture"
        ),
        **params,
    )

    model.load_state_dict(
        checkpoint["model_state_dict"],
        strict=True,
    )

    model.to(device)
    model.eval()

    return model

def _apply_label_inverse_transform(
    predictions: np.ndarray,
    task_names: list[str],
    label_scalers: Any,
    transform_metadata: Any,
) -> np.ndarray:
    """
    Convert model predictions from standardized transformed
    space back to original assay units.

    Training:
        raw value
        -> endpoint transform
        -> StandardScaler

    Inference:
        model output
        -> inverse StandardScaler
        -> inverse endpoint transform
    """
    predictions = np.asarray(
        predictions,
        dtype=np.float64,
    )

    if predictions.ndim != 2:
        raise ValueError(
            "Predictions must have shape "
            "[number_of_compounds, number_of_tasks]. "
            f"Received {predictions.shape}."
        )

    if predictions.shape[1] != len(
        task_names
    ):
        raise ValueError(
            "Prediction task dimension does not match "
            "the checkpoint task list. "
            f"Received {predictions.shape[1]} columns "
            f"for {len(task_names)} tasks."
        )

    restored = np.full(
        predictions.shape,
        np.nan,
        dtype=np.float64,
    )

    for task_index, task_name in enumerate(
        task_names
    ):
        if task_name not in label_scalers:
            raise KeyError(
                f"Checkpoint has no label scaler for "
                f"task '{task_name}'."
            )

        if task_name not in transform_metadata:
            raise KeyError(
                f"Checkpoint has no transform metadata "
                f"for task '{task_name}'."
            )

        scaler = label_scalers[
            task_name
        ]

        metadata = transform_metadata[
            task_name
        ]

        transform_name = metadata.get(
            "transform"
        )

        if transform_name is None:
            raise KeyError(
                f"Transform metadata for '{task_name}' "
                "does not contain 'transform'."
            )

        scaled_values = predictions[
            :,
            task_index,
        ]

        finite_mask = np.isfinite(
            scaled_values
        )

        if not finite_mask.any():
            continue

        transformed_values = np.full(
            scaled_values.shape,
            np.nan,
            dtype=np.float64,
        )

        transformed_values[
            finite_mask
        ] = (
            scaler
            .inverse_transform(
                scaled_values[
                    finite_mask
                ].reshape(-1, 1)
            )
            .reshape(-1)
        )

        if transform_name == "identity":
            assay_values = transformed_values

        elif transform_name == "log":
            assay_values = np.exp(
                transformed_values
            )

        elif transform_name == "log1p":
            assay_values = np.expm1(
                transformed_values
            )

        elif transform_name == "ic50_to_pic50":
            assay_values = pic50_to_ic50(
                transformed_values
            )

        elif transform_name == "ppb_to_logfu":
            assay_values = logfu_to_ppb(
                transformed_values
            )

        elif transform_name == "log_vd":
            assay_values = log_to_vd(
                transformed_values
            )

        elif transform_name == "logit_f":
            assay_values = (
                logit_to_bioavailability(
                    transformed_values
                )
            )

        else:
            raise ValueError(
                f"Unsupported transform "
                f"'{transform_name}' for task "
                f"'{task_name}'."
            )

        restored[
            :,
            task_index,
        ] = np.asarray(
            assay_values,
            dtype=np.float64,
        )

    return restored

class ADMECheckpointPredictor:
    """
    Framework-independent inference interface used by the PyFunc.
    """

    def __init__(
        self,
        checkpoint_path: str | Path,
        device: str | torch.device = "cpu",
        batch_size: int = 128,
    ):
        self.device = torch.device(
            device
        )

        self.batch_size = int(
            batch_size
        )

        if self.batch_size < 1:
            raise ValueError(
                "batch_size must be at least 1."
            )

        self.checkpoint = load_checkpoint(
            checkpoint_path
        )

        self.model = build_model_from_checkpoint(
            checkpoint=self.checkpoint,
            device=self.device,
        )

        self.task_names = list(
            self.checkpoint["task_names"]
        )

        self.label_scalers = self.checkpoint[
            "label_scalers"
        ]

        self.transform_metadata = self.checkpoint[
            "label_transform_metadata"
        ]

        self.global_feature_scaler = self.checkpoint[
            "global_feature_scaler"
        ]

    def _build_graphs(
        self,
        request_df: pd.DataFrame,
    ) -> tuple[
        list,
        list[int],
        list[str | None],
        list[str | None],
    ]:
        feature_state = prepare_features(
            df=request_df,
            smiles_col="smiles",
            global_scaler=(
                self.global_feature_scaler
            ),
            fit_global_scaler=False,
        )

        global_features = feature_state[
            "global_features"
        ]

        fp_features = feature_state[
            "fp_features"
        ]

        expected_rows = len(
            request_df
        )

        if len(global_features) != expected_rows:
            raise RuntimeError(
                "Global feature count does not match request size: "
                f"{len(global_features)} versus {expected_rows}."
            )

        if len(fp_features) != expected_rows:
            raise RuntimeError(
                "Fingerprint count does not match request size: "
                f"{len(fp_features)} versus {expected_rows}."
            )

        graphs = []
        valid_indices = []
        canonical_smiles = [
            None
        ] * expected_rows
        errors = [
            None
        ] * expected_rows

        # Use a local import to keep RDKit handling close to inference.
        from rdkit import Chem

        for position, smiles in enumerate(
            request_df["smiles"].tolist()
        ):

            if not smiles:
                errors[position] = (
                "SMILES is missing or empty."
                )
                continue
            mol = Chem.MolFromSmiles(
                smiles
            )

            if mol is None:
                errors[position] = (
                    "RDKit could not parse SMILES."
                )
                continue

            canonical_smiles[position] = (
                Chem.MolToSmiles(
                    mol,
                    canonical=True,
                    isomericSmiles=True,
                )
            )

            graph = mol_to_graph(
                smiles=smiles,
                label=None,
                idx=position,
                global_feats=(
                    global_features[position]
                ),
                fps=fp_features[position],
            )

            if graph is None:
                errors[position] = (
                    "Molecular graph construction failed."
                )
                continue

            graph.smiles = smiles

            graphs.append(
                graph
            )

            valid_indices.append(
                position
            )

        return (
            graphs,
            valid_indices,
            canonical_smiles,
            errors,
        )

    def predict(
        self,
        model_input: pd.DataFrame,
    ) -> pd.DataFrame:
        if not isinstance(
            model_input,
            pd.DataFrame,
        ):
            raise TypeError(
                "ADME prediction requires a pandas DataFrame."
            )

        if "smiles" not in model_input.columns:
            raise ValueError(
                "Input DataFrame requires a 'smiles' column."
            )

        request_df = (
            model_input
            .copy()
            .reset_index(drop=True)
        )

        if "compound_id" not in request_df.columns:
            request_df["compound_id"] = [
                f"row_{index}"
                for index in range(
                    len(request_df)
                )
            ]

        request_df["compound_id"] = (
            request_df["compound_id"]
            .astype(str)
        )

        missing_smiles = (
            request_df["smiles"].isna()
        )

        request_df.loc[
            missing_smiles,
            "smiles",
        ] = ""

        request_df["smiles"] = (
            request_df["smiles"]
            .astype(str)
            .str.strip()
        )

        output_df = pd.DataFrame(
            {
                "compound_id": request_df[
                    "compound_id"
                ],
                "input_smiles": request_df[
                    "smiles"
                ],
                "canonical_smiles": None,
                "valid_smiles": False,
                "error": None,
            }
        )

        for task_name in self.task_names:
            output_df[task_name] = np.nan

        if request_df.empty:
            return output_df

        (
            graphs,
            valid_indices,
            canonical_smiles,
            errors,
        ) = self._build_graphs(
            request_df
        )

        output_df["canonical_smiles"] = (
            canonical_smiles
        )

        output_df["error"] = errors

        if not graphs:
            return output_df

        loader = DataLoader(
            graphs,
            batch_size=self.batch_size,
            shuffle=False,
        )

        raw_batches = []

        with torch.inference_mode():
            for batch in loader:
                batch = batch.to(
                    self.device
                )

                raw_batch = self.model(
                    batch
                )

                raw_batches.append(
                    raw_batch
                    .detach()
                    .cpu()
                    .numpy()
                )

        raw_predictions = np.concatenate(
            raw_batches,
            axis=0,
        )

        expected_shape = (
            len(valid_indices),
            len(self.task_names),
        )

        if raw_predictions.shape != expected_shape:
            raise RuntimeError(
                "Unexpected raw prediction shape. "
                f"Expected {expected_shape}, "
                f"received {raw_predictions.shape}."
            )

        transformed_predictions = (
            _apply_label_inverse_transform(
                predictions=raw_predictions,
                task_names=self.task_names,
                label_scalers=(
                    self.label_scalers
                ),
                transform_metadata=(
                    self.transform_metadata
                ),
            )
        )

        transformed_predictions = np.asarray(
            transformed_predictions,
            dtype=np.float64,
        )

        if transformed_predictions.shape != expected_shape:
            raise RuntimeError(
                "Unexpected inverse-transformed prediction shape. "
                f"Expected {expected_shape}, received "
                f"{transformed_predictions.shape}."
            )

        for prediction_index, row_index in enumerate(
            valid_indices
        ):
            output_df.loc[
                row_index,
                "valid_smiles",
            ] = True

            output_df.loc[
                row_index,
                self.task_names,
            ] = transformed_predictions[
                prediction_index
            ]

        return output_df