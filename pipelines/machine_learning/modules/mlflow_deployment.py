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

# Root directory containing:
#   models/
#   modules/
#   pipeline/
#
# This file lives at:
#   pipelines/machine_learning/modules/mlflow_deployment.py
MACHINE_LEARNING_ROOT = (
    Path(__file__)
    .resolve()
    .parent
    .parent
)


def _resolve_code_paths(
    code_paths,
):
    """
    Resolve MLflow code paths relative to the
    machine_learning package root rather than
    the process working directory.
    """
    resolved_paths = []

    for code_path in code_paths:
        path = Path(
            os.path.expandvars(
                str(code_path)
            )
        ).expanduser()

        if not path.is_absolute():
            path = (
                MACHINE_LEARNING_ROOT
                / path
            )

        resolved_paths.append(
            str(path.resolve())
        )

    return resolved_paths

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

def _list_artifacts_recursive(
    client,
    run_id,
    path=None,
):
    results = []

    artifacts = client.list_artifacts(
        run_id=run_id,
        path=path,
    )

    for artifact in artifacts:
        results.append(
            artifact.path
        )

        if artifact.is_dir:
            results.extend(
                _list_artifacts_recursive(
                    client=client,
                    run_id=run_id,
                    path=artifact.path,
                )
            )

    return results

def _validate_code_paths(
    code_paths,
):
    """
    Verify that configured MLflow code paths exist
    and contain readable files before log_model()
    attempts to copy them.
    """
    logger.info(
        "Validating MLflow code paths..."
    )

    for code_path in code_paths:
        path = Path(
            code_path
        ).expanduser().resolve()

        logger.info(
            f"Checking code path: {path}"
        )

        if not path.exists():
            raise RuntimeError(
                "Configured MLflow code path "
                f"does not exist: {path}"
            )

        if path.is_file():
            try:
                with path.open("rb") as handle:
                    handle.read(1)
            except Exception as exc:
                raise RuntimeError(
                    "MLflow code file cannot be read: "
                    f"{path}: {exc}"
                ) from exc

            logger.info(
                f"Code file readable: {path}"
            )

            continue

        file_count = 0

        for item in path.rglob("*"):
            if not item.is_file():
                continue

            file_count += 1

            try:
                with item.open(
                    "rb"
                ) as handle:
                    handle.read(1)

            except Exception as exc:
                raise RuntimeError(
                    "MLflow code path contains "
                    "an unreadable file: "
                    f"{item}: {exc}"
                ) from exc

        logger.info(
            f"Code path readable: "
            f"{path} "
            f"({file_count} files)"
        )

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

    code_paths = _resolve_code_paths(
        pyfunc_cfg[
            "code_paths"
        ]
    )

    logger.info(
        f"Current working directory: "
        f"{Path.cwd()}"
    )

    logger.info(
        f"Machine-learning root: "
        f"{MACHINE_LEARNING_ROOT}"
    )

    logger.info(
        "Resolved MLflow code paths:"
    )

    for code_path in code_paths:
        logger.info(
            f"  {code_path}"
        )

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

        logger.info(
            f"Found migrated source run: "
            f"{source_run_id}"
        )

        source_run_info = mlflow.get_run(
            source_run_id
        )

        logger.info(
            f"Source artifact URI: "
            f"{source_run_info.info.artifact_uri}"
        )

        checkpoint_artifact_path = run_cfg[
            "checkpoint_artifact_path"
        ]

        logger.info(
            "Configured checkpoint artifact: "
            f"{checkpoint_artifact_path}"
        )

        # ---------------------------------
        # Inspect artifacts before download
        # ---------------------------------

        client = mlflow.MlflowClient()

        logger.info(
            "Inspecting source run artifacts..."
        )

        artifact_paths = (
            _list_artifacts_recursive(
                client=client,
                run_id=source_run_id,
            )
        )

        logger.info(
            f"Found {len(artifact_paths)} "
            "artifact path(s)"
        )

        for artifact_path in artifact_paths:
            logger.info(
                f"  {artifact_path}"
            )

        if (
            checkpoint_artifact_path
            not in artifact_paths
        ):
            raise RuntimeError(
                "Configured checkpoint artifact "
                "was not found.\n"
                f"Run ID: {source_run_id}\n"
                "Requested artifact: "
                f"{checkpoint_artifact_path}\n"
                "Available artifacts: "
                f"{artifact_paths}"
            )

        # ---------------------------------
        # Download verified artifact
        # ---------------------------------

        logger.info(
            "Checkpoint artifact exists. "
            "Starting download..."
        )

        checkpoint_path = (
            mlflow.artifacts.download_artifacts(
                run_id=source_run_id,
                artifact_path=(
                    checkpoint_artifact_path
                ),
            )
        )

        logger.info(
            "Checkpoint download complete."
        )

        checkpoint_file = Path(
            checkpoint_path
        )

        if not checkpoint_file.is_file():
            raise RuntimeError(
                "Downloaded checkpoint is not "
                f"a file: {checkpoint_file}"
            )

        logger.info(
            "Checkpoint verified locally: "
            f"{checkpoint_file}"
        )

        logger.info(
            "Beginning PyFunc packaging..."
        )

        # checkpoint_candidate = Path(checkpoint_path)

        # if checkpoint_candidate.is_file():
        #     checkpoint_file = checkpoint_candidate
        # elif checkpoint_candidate.is_dir():
        #     checkpoint_files = sorted(
        #         checkpoint_candidate.rglob("*.pth")
        #     )

        #     if not checkpoint_files:
        #         raise RuntimeError(
        #             f"No checkpoint found for source run "
        #             f"{source_run_id} under {checkpoint_candidate}"
        #         )

        #     if len(checkpoint_files) > 1:
        #         raise RuntimeError(
        #             "Multiple checkpoints found. Configure the exact "
        #             f"artifact path instead: {checkpoint_files}"
        #         )

        #     checkpoint_file = checkpoint_files[0]
        # else:
        #     raise RuntimeError(
        #         f"Checkpoint artifact does not exist: "
        #         f"{checkpoint_candidate}"
        #     )

        # checkpoint_file = Path(
        #     checkpoint_path
        # )

        # if not checkpoint_file.is_file():
        #     raise RuntimeError(
        #         "Downloaded checkpoint artifact is not "
        #         f"a file: {checkpoint_file}"
        #     )

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

            _validate_code_paths(
                code_paths
            )

            logger.info(
                "Code paths validated. "
                "Calling mlflow.pyfunc.log_model..."
            )
        
            try:
                model_info = mlflow.pyfunc.log_model(
                    name=pyfunc_cfg[
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
                        "checkpoint":
                            checkpoint_path,
                    },

                    code_paths=code_paths,

                    pip_requirements=pyfunc_cfg[
                        "pip_requirements"
                    ],

                    input_example=input_example,
                )

            except Exception:
                logger.exception(
                    "MLflow PyFunc logging failed."
                )
                raise

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
