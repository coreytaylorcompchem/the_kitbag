from __future__ import annotations

import mlflow.pyfunc
import pandas as pd
import torch

from models.adme.mlflow_inference import (
    ADMECheckpointPredictor,
)


class ADMEMultitaskPyFunc(
    mlflow.pyfunc.PythonModel
):
    """
    MLflow PythonModel wrapper for multitask ADME inference.
    """

    def __init__(
        self,
        batch_size: int = 128,
        prefer_gpu: bool = False,
    ):
        self.batch_size = int(
            batch_size
        )

        self.prefer_gpu = bool(
            prefer_gpu
        )

        self.predictor = None

    def load_context(
        self,
        context,
    ) -> None:
        use_gpu = (
            self.prefer_gpu
            and torch.cuda.is_available()
        )

        device = (
            "cuda"
            if use_gpu
            else "cpu"
        )

        self.predictor = ADMECheckpointPredictor(
            checkpoint_path=(
                context.artifacts[
                    "checkpoint"
                ]
            ),
            device=device,
            batch_size=self.batch_size,
        )

    def predict(
        self,
        context,
        model_input: pd.DataFrame,
        params=None,
    ) -> pd.DataFrame:
        if self.predictor is None:
            raise RuntimeError(
                "ADME predictor has not been initialised."
            )

        return self.predictor.predict(
            model_input
        )