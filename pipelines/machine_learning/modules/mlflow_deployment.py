from pathlib import Path
import os

import mlflow
import pandas as pd

from pipeline.logger import setup_logger
from pipeline.task_registry import register_task

from models.adme.mtl_adme.mlflow_pyfunc import (
ADMEMultitaskPyFunc,
)

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True,
)

def _find_migrated_run(
    experiment_name,
    model_id,
):
    experiment = mlflow.get_experiment_by_name(
        experiment_name
    )

    if experiment is None:
        raise RuntimeError(
            f"Experiment not found: {experiment_name}"
        )

    runs = mlflow.search_runs(
        experiment_ids=[
            experiment.experiment_id
        ],
        filter_string=(
            f"tags.`migration.model_id` = '{model_id}'"
        ),
    )

    if len(runs) == 0:
        raise RuntimeError(
            f"No migrated run found for "
            f"model_id={model_id}"
        )

    if len(runs) > 1:
        raise RuntimeError(
            f"Found {len(runs)} migrated runs for "
            f"model_id={model_id}. "
            "Expected exactly one source run."
        )

    return runs.iloc[0]

@register_task(
    "create_mlflow_serving_models",
    category="MLflow",
)
def create_mlflow_serving_models(
    config,
    context,
):

    tracking_uri = os.path.expandvars(
        config["tracking_uri"]
    )

    mlflow.set_tracking_uri(
        tracking_uri
    )

    experiment_name = config[
        "experiment_name"
    ]

    run_configs = config[
        "runs"
    ]

    pyfunc_cfg = config[
        "pyfunc"
    ]

    deployment_results = []

    for run_cfg in run_configs:

        model_id = run_cfg[
            "model_id"
        ]

        logger.info(
            f"Deploying {model_id}"
        )

        source_run = _find_migrated_run(
            experiment_name,
            model_id,
        )

        source_run_id = source_run[
            "run_id"
        ]

        checkpoint_path = mlflow.artifacts.download_artifacts(
            run_id=source_run_id,
            artifact_path=run_cfg[
                "checkpoint_artifact_path"
            ],
        )

        checkpoint_candidate = Path(checkpoint_path)

        if checkpoint_candidate.is_file():
            checkpoint_file = checkpoint_candidate
        elif checkpoint_candidate.is_dir():
            checkpoint_files = sorted(
                checkpoint_candidate.rglob("*.pth")
            )

            if not checkpoint_files:
                raise RuntimeError(
                    f"No checkpoint found for source run "
                    f"{source_run_id} under {checkpoint_candidate}"
                )

            if len(checkpoint_files) > 1:
                raise RuntimeError(
                    "Multiple checkpoints found. Configure the exact "
                    f"artifact path instead: {checkpoint_files}"
                )

            checkpoint_file = checkpoint_files[0]
        else:
            raise RuntimeError(
                f"Checkpoint artifact does not exist: "
                f"{checkpoint_candidate}"
            )

        checkpoint_file = Path(
            checkpoint_path
        )

        if not checkpoint_file.is_file():
            raise RuntimeError(
                "Downloaded checkpoint artifact is not "
                f"a file: {checkpoint_file}"
            )

        checkpoint_path = str(
            checkpoint_file
        )

        experiment = mlflow.get_experiment_by_name(
            experiment_name
        )

        with mlflow.start_run(
            experiment_id=experiment.experiment_id,
            run_name=f"deploy_{model_id}",
        ):

            mlflow.set_tags(
                {
                    "deployment": "pyfunc",
                    "model_id": model_id,
                    "source_run_id": source_run_id,
                }
            )

            input_example = pd.DataFrame(
                {
                    "compound_id": [
                        "example"
                    ],
                    "smiles": [
                        "CCO"
                    ],
                }
            )
        
            model_info = mlflow.pyfunc.log_model(
                artifact_path=pyfunc_cfg[
                    "model_name"
                ],

                python_model=ADMEMultitaskPyFunc(
                    batch_size=pyfunc_cfg.get(
                        "batch_size",
                        128,
                    ),
                    prefer_gpu=pyfunc_cfg.get(
                        "prefer_gpu",
                        False,
                    ),
                ),

                artifacts={
                    "checkpoint": checkpoint_path,
                },

                code_paths=pyfunc_cfg[
                    "code_paths"
                ],

                pip_requirements=pyfunc_cfg[
                    "pip_requirements"
                ],

                input_example=input_example,
            )
        
            deployment_results.append(
                {
                    "model_id": model_id,
                    "source_run_id": source_run_id,
                    "deployment_run_id":
                        mlflow.active_run().info.run_id,
                    "model_uri":
                        model_info.model_uri,
                }
            )

    logger.info(
        f"Created "
        f"{len(deployment_results)} "
        f"serving model(s)"
    )

    return {
        "deployment_results":
            deployment_results
    }
