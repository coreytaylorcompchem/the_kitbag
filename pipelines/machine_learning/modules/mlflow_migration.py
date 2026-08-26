import hashlib
import json
import math
import os
import platform
import re
import socket
import sys

from dotenv import load_dotenv
from pathlib import Path
from typing import Any, Optional

import mlflow
import pandas as pd
import torch
import yaml

from mlflow.tracking import MlflowClient
from mlflow.models import ModelSignature
from mlflow.types import ColSpec, Schema

# from inference.adme_pyfunc_model import ADMEMultitaskPyFunc

from pipeline.logger import setup_logger
from pipeline.task_registry import register_task

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True,
)


# =============================================================================
# General helpers
# =============================================================================

def _validate_mlflow_authentication_environment() -> None:
    """
    Validate that MLflow username/password credentials are available.

    Credentials are intentionally never logged.
    """
    username = os.environ.get("MLFLOW_TRACKING_USERNAME")
    password = os.environ.get("MLFLOW_TRACKING_PASSWORD")

    missing = []

    if not username:
        missing.append("MLFLOW_TRACKING_USERNAME")

    if not password:
        missing.append("MLFLOW_TRACKING_PASSWORD")

    if missing:
        raise RuntimeError(
            "MLflow authentication credentials are missing. "
            f"Expected environment variables: {', '.join(missing)}. "
            "These can be supplied through a local .env file."
        )

    logger.info(
        "MLflow username/password authentication credentials detected."
    )

def _resolve_environment_value(value: Any) -> Any:
    """
    Resolve YAML values written as:

        ${ENVIRONMENT_VARIABLE}

    This is included in case the pipeline YAML loader does not already perform
    environment-variable interpolation.
    """
    if not isinstance(value, str):
        return value

    match = re.fullmatch(
        r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}",
        value.strip(),
    )

    if match is None:
        return value

    variable_name = match.group(1)
    resolved_value = os.environ.get(variable_name)

    if resolved_value is None:
        raise RuntimeError(
            f"Environment variable '{variable_name}' is required but is not set."
        )

    return resolved_value


def _as_path(value: Any) -> Optional[Path]:
    if value is None:
        return None

    if isinstance(value, float) and pd.isna(value):
        return None

    text = str(value).strip()

    if not text:
        return None

    return Path(text).expanduser().resolve()


def _sanitize_mlflow_key(
    value: Any,
    max_length: int = 240,
) -> str:
    """
    Convert task and metric names to conservative MLflow-compatible keys.
    """
    text = str(value).strip()

    text = re.sub(r"\s+", "_", text)
    text = re.sub(r"[^A-Za-z0-9_.\-/]", "_", text)
    text = re.sub(r"_+", "_", text)
    text = text.strip("._-/")

    if not text:
        text = "unnamed"

    return text[:max_length]


def _sanitize_run_name(value: Any) -> str:
    text = str(value).strip()

    text = re.sub(r"\s+", "_", text)
    text = re.sub(r"[^A-Za-z0-9_.-]", "_", text)
    text = re.sub(r"_+", "_", text)

    return text.strip("._-") or "historical_model"


def _sha256_file(
    path: Path,
    chunk_size: int = 8 * 1024 * 1024,
) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)

            if not chunk:
                break

            digest.update(chunk)

    return digest.hexdigest()


def _json_safe(value: Any) -> Any:
    if value is None:
        return None

    if isinstance(value, Path):
        return str(value)

    if isinstance(value, dict):
        return {
            str(key): _json_safe(item)
            for key, item in value.items()
        }

    if isinstance(value, (list, tuple, set)):
        return [
            _json_safe(item)
            for item in value
        ]

    if isinstance(value, (str, int, float, bool)):
        return value

    return str(value)


def _load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)

    if data is None:
        return {}

    if not isinstance(data, dict):
        raise TypeError(
            f"Expected YAML root to be a mapping: {path}"
        )

    return data


def _flatten_mapping(
    mapping: dict[str, Any],
    prefix: str = "",
) -> dict[str, Any]:
    """
    Flatten nested YAML values for MLflow parameter logging.

    Lists are stored as compact JSON strings. The complete source YAML is also
    uploaded as an artifact, so MLflow parameters are only a searchable view.
    """
    flattened: dict[str, Any] = {}

    for key, value in mapping.items():
        full_key = f"{prefix}.{key}" if prefix else str(key)
        full_key = _sanitize_mlflow_key(full_key)

        if isinstance(value, dict):
            flattened.update(
                _flatten_mapping(
                    value,
                    prefix=full_key,
                )
            )

        elif isinstance(value, (list, tuple)):
            flattened[full_key] = json.dumps(
                _json_safe(value),
                separators=(",", ":"),
            )

        elif value is None:
            flattened[full_key] = "null"

        else:
            flattened[full_key] = value

    return flattened


# =============================================================================
# MLflow logging helpers
# =============================================================================

def _safe_log_params(
    params: dict[str, Any],
    batch_size: int = 100,
) -> None:
    """
    Log parameters in batches.

    The complete configuration remains available as an uploaded YAML artifact.
    Long parameter values are truncated only in MLflow's searchable parameter
    view.
    """
    cleaned: dict[str, str] = {}

    for key, value in params.items():
        clean_key = _sanitize_mlflow_key(key)
        clean_value = str(_json_safe(value))

        if len(clean_value) > 500:
            clean_value = clean_value[:497] + "..."

        cleaned[clean_key] = clean_value

    items = list(cleaned.items())

    for start in range(0, len(items), batch_size):
        batch = dict(items[start:start + batch_size])

        if batch:
            mlflow.log_params(batch)


def _safe_log_metrics(
    metrics: dict[str, float],
    batch_size: int = 100,
) -> None:
    """
    Log finite numeric metrics in batches.
    """
    cleaned: dict[str, float] = {}

    for key, value in metrics.items():
        try:
            numeric_value = float(value)
        except (TypeError, ValueError):
            continue

        if not math.isfinite(numeric_value):
            continue

        cleaned[_sanitize_mlflow_key(key)] = numeric_value

    items = list(cleaned.items())

    for start in range(0, len(items), batch_size):
        batch = dict(items[start:start + batch_size])

        if batch:
            mlflow.log_metrics(batch)


def _runtime_metadata() -> dict[str, Any]:
    return {
        "python_version": sys.version,
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "hostname": socket.gethostname(),
        "torch_version": torch.__version__,
        "cuda_available": torch.cuda.is_available(),
        "cuda_version": torch.version.cuda,
        "mlflow_version": mlflow.__version__,
        "pandas_version": pd.__version__,
    }


# =============================================================================
# Metrics CSV handling
# =============================================================================

def _find_task_column(
    metrics_df: pd.DataFrame,
    candidates: list[str],
) -> str:
    lower_to_original = {
        str(column).strip().lower(): str(column)
        for column in metrics_df.columns
    }

    for candidate in candidates:
        candidate_lower = str(candidate).strip().lower()

        if candidate_lower in lower_to_original:
            return lower_to_original[candidate_lower]

    non_numeric_columns = [
        str(column)
        for column in metrics_df.columns
        if not pd.api.types.is_numeric_dtype(metrics_df[column])
    ]

    if len(non_numeric_columns) == 1:
        inferred_column = non_numeric_columns[0]

        logger.info(
            f"Inferred metrics task column: {inferred_column}"
        )

        return inferred_column

    raise ValueError(
        "Could not identify the task column in the metrics CSV. "
        f"Candidate names were: {candidates}. "
        f"Available columns are: {list(metrics_df.columns)}"
    )


def _find_metric_columns(
    metrics_df: pd.DataFrame,
    task_column: str,
    configured_columns: list[str],
) -> list[str]:
    if configured_columns:
        missing_columns = [
            column
            for column in configured_columns
            if column not in metrics_df.columns
        ]

        if missing_columns:
            raise KeyError(
                "Configured metrics columns are missing from the metrics CSV: "
                f"{missing_columns}. "
                f"Available columns: {list(metrics_df.columns)}"
            )

        return list(configured_columns)

    numeric_columns: list[str] = []

    for column in metrics_df.columns:
        if column == task_column:
            continue

        converted = pd.to_numeric(
            metrics_df[column],
            errors="coerce",
        )

        if converted.notna().any():
            numeric_columns.append(str(column))

    if not numeric_columns:
        raise ValueError(
            "No numeric metric columns were found in the metrics CSV. "
            f"Available columns: {list(metrics_df.columns)}"
        )

    return numeric_columns


def _extract_mlflow_metrics(
    metrics_df: pd.DataFrame,
    metrics_config: dict[str, Any],
) -> tuple[
    str,
    list[str],
    dict[str, float],
    dict[str, float],
]:
    candidates = metrics_config.get(
        "task_column_candidates",
        [
            "task",
            "task_name",
            "endpoint",
            "property",
            "target",
            "label",
        ],
    )

    task_column = _find_task_column(
        metrics_df=metrics_df,
        candidates=candidates,
    )

    configured_metric_columns = metrics_config.get(
        "metric_columns",
        [],
    )

    metric_columns = _find_metric_columns(
        metrics_df=metrics_df,
        task_column=task_column,
        configured_columns=configured_metric_columns,
    )

    metric_prefix = _sanitize_mlflow_key(
        metrics_config.get("metric_prefix", "task")
    )

    summary_prefix = _sanitize_mlflow_key(
        metrics_config.get("summary_prefix", "summary")
    )

    per_task_metrics: dict[str, float] = {}

    duplicate_task_names = (
        metrics_df[task_column]
        .astype(str)
        .duplicated()
        .sum()
    )

    if duplicate_task_names > 0:
        raise ValueError(
            f"Metrics CSV contains {duplicate_task_names} duplicate task rows "
            f"in column '{task_column}'. Expected one row per task."
        )

    for _, row in metrics_df.iterrows():
        raw_task_name = row[task_column]

        if pd.isna(raw_task_name):
            logger.warning(
                "Skipping metrics row with missing task name."
            )
            continue

        task_name = _sanitize_mlflow_key(raw_task_name)

        for metric_column in metric_columns:
            numeric_value = pd.to_numeric(
                pd.Series([row[metric_column]]),
                errors="coerce",
            ).iloc[0]

            if pd.isna(numeric_value):
                continue

            numeric_value = float(numeric_value)

            if not math.isfinite(numeric_value):
                continue

            metric_key = (
                f"{metric_prefix}."
                f"{task_name}."
                f"{_sanitize_mlflow_key(metric_column)}"
            )

            per_task_metrics[metric_key] = numeric_value

    summary_metrics: dict[str, float] = {
        f"{summary_prefix}.n_tasks": float(
            metrics_df[task_column].nunique()
        )
    }

    for metric_column in metric_columns:
        values = pd.to_numeric(
            metrics_df[metric_column],
            errors="coerce",
        )

        values = values[
            values.notna()
            & values.map(
                lambda value: math.isfinite(float(value))
            )
        ]

        if values.empty:
            continue

        metric_name = _sanitize_mlflow_key(metric_column)

        summary_metrics[
            f"{summary_prefix}.macro_mean.{metric_name}"
        ] = float(values.mean())

        summary_metrics[
            f"{summary_prefix}.macro_median.{metric_name}"
        ] = float(values.median())

        summary_metrics[
            f"{summary_prefix}.n_valid.{metric_name}"
        ] = float(len(values))

    return (
        task_column,
        metric_columns,
        per_task_metrics,
        summary_metrics,
    )


# =============================================================================
# Historical model validation and artifact logging
# =============================================================================

def _validate_model_config(
    model_config: dict[str, Any],
    require_optional_artifacts: bool,
) -> dict[str, Any]:
    required_keys = [
        "model_id",
        "checkpoint_path",
        "metrics_csv_path",
    ]

    missing_keys = [
        key
        for key in required_keys
        if not model_config.get(key)
    ]

    if missing_keys:
        raise KeyError(
            f"Historical model block is missing required keys: {missing_keys}"
        )

    resolved = dict(model_config)

    path_keys = [
        "checkpoint_path",
        "metrics_csv_path",
        "scaler_path",
        "transform_metadata_path",
        "training_config_path",
    ]

    for key in path_keys:
        resolved[key] = _as_path(model_config.get(key))

    resolved["artifact_dirs"] = [
        _as_path(path)
        for path in model_config.get("artifact_dirs", [])
        if path is not None
    ]

    checkpoint_path = resolved["checkpoint_path"]
    metrics_csv_path = resolved["metrics_csv_path"]

    if checkpoint_path is None or not checkpoint_path.is_file():
        raise FileNotFoundError(
            f"Checkpoint does not exist: {checkpoint_path}"
        )

    if metrics_csv_path is None or not metrics_csv_path.is_file():
        raise FileNotFoundError(
            f"Metrics CSV does not exist: {metrics_csv_path}"
        )

    optional_file_keys = [
        "scaler_path",
        "transform_metadata_path",
        "training_config_path",
    ]

    for key in optional_file_keys:
        path = resolved.get(key)

        if path is None:
            continue

        if not path.is_file():
            message = (
                f"Optional artifact configured for model "
                f"'{resolved['model_id']}' does not exist: {path}"
            )

            if require_optional_artifacts:
                raise FileNotFoundError(message)

            logger.warning(message)
            resolved[key] = None

    valid_artifact_dirs = []

    for artifact_dir in resolved["artifact_dirs"]:
        if artifact_dir is None:
            continue

        if not artifact_dir.is_dir():
            message = (
                f"Artifact directory configured for model "
                f"'{resolved['model_id']}' does not exist: {artifact_dir}"
            )

            if require_optional_artifacts:
                raise NotADirectoryError(message)

            logger.warning(message)
            continue

        valid_artifact_dirs.append(artifact_dir)

    resolved["artifact_dirs"] = valid_artifact_dirs

    return resolved


def _log_model_artifacts(
    model_config: dict[str, Any],
) -> None:
    checkpoint_path = model_config["checkpoint_path"]
    metrics_csv_path = model_config["metrics_csv_path"]

    logger.info(f"Uploading checkpoint: {checkpoint_path}")

    mlflow.log_artifact(
        str(checkpoint_path),
        artifact_path="model/checkpoint",
    )

    logger.info(f"Uploading metrics CSV: {metrics_csv_path}")

    mlflow.log_artifact(
        str(metrics_csv_path),
        artifact_path="evaluation",
    )

    optional_files = [
        (
            "scaler_path",
            "model/preprocessing",
        ),
        (
            "transform_metadata_path",
            "model/preprocessing",
        ),
        (
            "training_config_path",
            "configuration",
        ),
    ]

    for config_key, artifact_path in optional_files:
        path = model_config.get(config_key)

        if path is None:
            continue

        logger.info(f"Uploading artifact: {path}")

        mlflow.log_artifact(
            str(path),
            artifact_path=artifact_path,
        )

    for index, artifact_dir in enumerate(
        model_config.get("artifact_dirs", []),
        start=1,
    ):
        directory_name = _sanitize_mlflow_key(
            artifact_dir.name
        )

        remote_artifact_path = (
            f"evaluation/artifact_directories/"
            f"{index:02d}_{directory_name}"
        )

        logger.info(
            f"Uploading artifact directory: {artifact_dir}"
        )

        mlflow.log_artifacts(
            str(artifact_dir),
            artifact_path=remote_artifact_path,
        )

def _log_deployable_model(
    model_config: dict[str, Any],
) -> None:
    checkpoint_path = model_config["checkpoint_path"]
    inference_config_path = model_config.get(
        "inference_config_path"
    )

    if inference_config_path is None:
        raise ValueError(
            "A deployable model requires inference_config_path."
        )

    artifacts = {
        "checkpoint": str(checkpoint_path),
        "inference_config": str(inference_config_path),
    }

    scaler_path = model_config.get("scaler_path")

    if scaler_path is not None:
        artifacts["label_scalers"] = str(scaler_path)

    transform_metadata_path = model_config.get(
        "transform_metadata_path"
    )

    if transform_metadata_path is not None:
        artifacts["transform_metadata"] = str(
            transform_metadata_path
        )

    input_schema = Schema(
        [
            ColSpec("string", "compound_id", required=False),
            ColSpec("string", "smiles", required=True),
        ]
    )

    output_schema = Schema(
        [
            ColSpec("string", "compound_id"),
            ColSpec("string", "input_smiles"),
            ColSpec("string", "canonical_smiles"),
            ColSpec("boolean", "valid_smiles"),
            ColSpec("string", "error"),
        ]
    )

    signature = ModelSignature(
        inputs=input_schema,
        outputs=output_schema,
    )

    input_example = pd.DataFrame(
        {
            "compound_id": ["example_1", "example_2"],
            "smiles": [
                "CCO",
                "CC(=O)Oc1ccccc1C(=O)O",
            ],
        }
    )

    mlflow.pyfunc.log_model(
        artifact_path="deployable_model",
        python_model=ADMEMultitaskPyFunc(),
        artifacts=artifacts,
        code_paths=[
            "inference/adme_pyfunc_model.py",
            "inference/adme_inference_adapter.py",
            "models",
        ],
        pip_requirements=[
            f"mlflow=={mlflow.__version__}",
            f"torch=={torch.__version__}",
            "pandas",
            "numpy",
            "rdkit",
            "torch-geometric",
            "joblib",
            "pyyaml",
        ],
        signature=signature,
        input_example=input_example,
    )

def _build_import_manifest(
    model_config: dict[str, Any],
    run_id: str,
    checkpoint_sha256: Optional[str],
    metrics_sha256: Optional[str],
    task_column: str,
    metric_columns: list[str],
    metrics_df: pd.DataFrame,
) -> dict[str, Any]:
    checkpoint_path = model_config["checkpoint_path"]
    metrics_csv_path = model_config["metrics_csv_path"]

    return {
        "model_id": model_config["model_id"],
        "run_id": run_id,
        "source": "historical_model_migration",
        "checkpoint": {
            "path": str(checkpoint_path),
            "filename": checkpoint_path.name,
            "size_bytes": checkpoint_path.stat().st_size,
            "sha256": checkpoint_sha256,
        },
        "metrics": {
            "path": str(metrics_csv_path),
            "filename": metrics_csv_path.name,
            "size_bytes": metrics_csv_path.stat().st_size,
            "sha256": metrics_sha256,
            "task_column": task_column,
            "metric_columns": metric_columns,
            "number_of_rows": int(len(metrics_df)),
            "number_of_tasks": int(
                metrics_df[task_column].nunique()
            ),
        },
        "additional_artifacts": {
            "scaler_path": (
                str(model_config["scaler_path"])
                if model_config.get("scaler_path") is not None
                else None
            ),
            "transform_metadata_path": (
                str(model_config["transform_metadata_path"])
                if model_config.get("transform_metadata_path") is not None
                else None
            ),
            "training_config_path": (
                str(model_config["training_config_path"])
                if model_config.get("training_config_path") is not None
                else None
            ),
            "artifact_dirs": [
                str(path)
                for path in model_config.get("artifact_dirs", [])
            ],
        },
        "runtime": _runtime_metadata(),
    }


# =============================================================================
# Idempotency
# =============================================================================

def _escape_mlflow_filter_value(value: str) -> str:
    return value.replace("\\", "\\\\").replace("'", "\\'")


def _find_existing_run(
    client: MlflowClient,
    experiment_id: str,
    model_id: str,
    checkpoint_sha256: str,
) -> Optional[str]:
    escaped_model_id = _escape_mlflow_filter_value(model_id)
    escaped_sha256 = _escape_mlflow_filter_value(checkpoint_sha256)

    filter_string = (
        f"tags.`migration.model_id` = '{escaped_model_id}' "
        f"AND "
        f"tags.`migration.checkpoint_sha256` = '{escaped_sha256}' "
        f"AND "
        f"tags.`migration.status` = 'completed'"
    )

    matching_runs = client.search_runs(
        experiment_ids=[experiment_id],
        filter_string=filter_string,
        max_results=1,
    )

    if not matching_runs:
        return None

    return matching_runs[0].info.run_id


# =============================================================================
# Single-model migration
# =============================================================================

def _migrate_single_model(
    model_config: dict[str, Any],
    experiment_id: str,
    metrics_config: dict[str, Any],
    skip_existing: bool,
    calculate_checksums: bool,
    require_optional_artifacts: bool,
) -> dict[str, Any]:

    model_config = _validate_model_config(
        model_config=model_config,
        require_optional_artifacts=require_optional_artifacts,
    )

    model_id = str(model_config["model_id"])
    checkpoint_path = model_config["checkpoint_path"]
    metrics_csv_path = model_config["metrics_csv_path"]

    checkpoint_sha256 = None
    metrics_sha256 = None

    if calculate_checksums:
        logger.info(
            f"Calculating checkpoint SHA256 for model '{model_id}'"
        )

        checkpoint_sha256 = _sha256_file(checkpoint_path)

        logger.info(
            f"Calculating metrics CSV SHA256 for model '{model_id}'"
        )

        metrics_sha256 = _sha256_file(metrics_csv_path)

    client = MlflowClient()

    if (
        skip_existing
        and checkpoint_sha256 is not None
    ):
        existing_run_id = _find_existing_run(
            client=client,
            experiment_id=experiment_id,
            model_id=model_id,
            checkpoint_sha256=checkpoint_sha256,
        )

        if existing_run_id is not None:
            logger.info(
                f"Skipping model '{model_id}'. "
                f"The checkpoint was already imported in run "
                f"'{existing_run_id}'."
            )

            return {
                "model_id": model_id,
                "status": "skipped_existing",
                "run_id": existing_run_id,
                "checkpoint_path": str(checkpoint_path),
                "metrics_csv_path": str(metrics_csv_path),
                "checkpoint_sha256": checkpoint_sha256,
                "error": None,
            }

    logger.info(f"Reading metrics CSV: {metrics_csv_path}")

    metrics_df = pd.read_csv(metrics_csv_path)

    (
        task_column,
        metric_columns,
        per_task_metrics,
        summary_metrics,
    ) = _extract_mlflow_metrics(
        metrics_df=metrics_df,
        metrics_config=metrics_config,
    )

    run_name = _sanitize_run_name(
        model_config.get(
            "run_name",
            f"historical_{model_id}",
        )
    )

    user_tags = {
        _sanitize_mlflow_key(key): str(value)
        for key, value in model_config.get("tags", {}).items()
    }

    run_tags = {
        "migration.model_id": model_id,
        "migration.source": "historical_checkpoint",
        "migration.status": "in_progress",
        "model.framework": "pytorch",
        **user_tags,
    }

    if checkpoint_sha256 is not None:
        run_tags["migration.checkpoint_sha256"] = checkpoint_sha256

    if metrics_sha256 is not None:
        run_tags["migration.metrics_sha256"] = metrics_sha256

    with mlflow.start_run(
        experiment_id=experiment_id,
        run_name=run_name,
        tags=run_tags,
    ) as active_run:

        run_id = active_run.info.run_id

        try:
            logger.info(
                f"Created MLflow run '{run_id}' for model '{model_id}'"
            )

            base_params = {
                "model_id": model_id,
                "historical_import": True,
                "checkpoint_filename": checkpoint_path.name,
                "metrics_filename": metrics_csv_path.name,
                "number_of_tasks": metrics_df[task_column].nunique(),
                "number_of_metric_columns": len(metric_columns),
            }

            _safe_log_params(base_params)

            training_config_path = model_config.get(
                "training_config_path"
            )

            if training_config_path is not None:
                logger.info(
                    f"Logging searchable parameters from "
                    f"{training_config_path}"
                )

                training_config = _load_yaml(
                    training_config_path
                )

                flattened_training_config = _flatten_mapping(
                    training_config
                )

                prefixed_training_params = {
                    f"training_config.{key}": value
                    for key, value in flattened_training_config.items()
                }

                _safe_log_params(
                    prefixed_training_params
                )

            logger.info(
                f"Logging {len(summary_metrics)} summary metrics"
            )

            _safe_log_metrics(summary_metrics)

            logger.info(
                f"Logging {len(per_task_metrics)} per-task metrics"
            )

            _safe_log_metrics(per_task_metrics)

            _log_model_artifacts(model_config)

            import_manifest = _build_import_manifest(
                model_config=model_config,
                run_id=run_id,
                checkpoint_sha256=checkpoint_sha256,
                metrics_sha256=metrics_sha256,
                task_column=task_column,
                metric_columns=metric_columns,
                metrics_df=metrics_df,
            )

            mlflow.log_dict(
                _json_safe(import_manifest),
                "migration/import_manifest.json",
            )

            mlflow.log_dict(
                _json_safe(_runtime_metadata()),
                "environment/runtime.json",
            )

            mlflow.set_tags(
                {
                    "migration.status": "completed",
                    "migration.n_tasks": str(
                        metrics_df[task_column].nunique()
                    ),
                    "migration.n_logged_metrics": str(
                        len(per_task_metrics)
                        + len(summary_metrics)
                    ),
                }
            )

            logger.info(
                f"Completed MLflow migration for model '{model_id}'"
            )

            return {
                "model_id": model_id,
                "status": "completed",
                "run_id": run_id,
                "checkpoint_path": str(checkpoint_path),
                "metrics_csv_path": str(metrics_csv_path),
                "checkpoint_sha256": checkpoint_sha256,
                "n_tasks": int(
                    metrics_df[task_column].nunique()
                ),
                "n_per_task_metrics": len(per_task_metrics),
                "n_summary_metrics": len(summary_metrics),
                "error": None,
            }

        except Exception as error:
            mlflow.set_tags(
                {
                    "migration.status": "failed",
                    "migration.error_type": type(error).__name__,
                    "migration.error": str(error)[:5000],
                }
            )

            raise


# =============================================================================
# Registered pipeline task
# =============================================================================

@register_task(
    "migrate_adme_models_to_mlflow",
    category="ADME",
    description=(
        "Migrate existing ADME checkpoints, metrics, configuration, "
        "and evaluation artifacts to MLflow without retraining."
    ),
)
def migrate_adme_models_to_mlflow(config, context):
    """
    Migrate historical ADME models to a remote MLflow tracking server.

    Expected YAML structure:

        migrate_adme_models_to_mlflow:
          tracking_uri: ${MLFLOW_TRACKING_URI}
          experiment_name: ADME_historical_models
          skip_existing: true
          continue_on_error: true
          calculate_checksums: true
          migration_report_path: outputs/mlflow/report.csv
          metrics: ...
          models:
            - model_id: ...
              checkpoint_path: ...
              metrics_csv_path: ...

    One MLflow run is created per historical model.
    """

    tracking_uri = _resolve_environment_value(
        config.get("tracking_uri")
    )

    if not tracking_uri:
        raise ValueError(
            "migrate_adme_models_to_mlflow requires 'tracking_uri'. "
            "Set it directly in YAML or use ${MLFLOW_TRACKING_URI}."
        )

    experiment_name = config.get(
        "experiment_name",
        "ADME_historical_models",
    )

    models = config.get("models", [])

    if not models:
        raise ValueError(
            "migrate_adme_models_to_mlflow requires a non-empty "
            "'models' list."
        )

    skip_existing = config.get(
        "skip_existing",
        True,
    )

    continue_on_error = config.get(
        "continue_on_error",
        True,
    )

    calculate_checksums = config.get(
        "calculate_checksums",
        True,
    )

    require_optional_artifacts = config.get(
        "require_optional_artifacts",
        False,
    )

    metrics_config = config.get(
        "metrics",
        {},
    )

    report_path = Path(
        config.get(
            "migration_report_path",
            "outputs/mlflow/adme_model_migration_report.csv",
        )
    ).expanduser().resolve()

    logger.info("=" * 80)
    logger.info("Starting historical ADME model migration")
    logger.info(f"MLflow tracking URI: {tracking_uri}")
    logger.info(f"MLflow experiment: {experiment_name}")
    logger.info(f"Number of configured models: {len(models)}")
    logger.info(f"Skip existing checkpoints: {skip_existing}")
    logger.info(f"Continue on error: {continue_on_error}")
    logger.info("=" * 80)

    _validate_mlflow_authentication_environment()

    mlflow.set_tracking_uri(tracking_uri)

    experiment = mlflow.set_experiment(
        experiment_name
    )

    experiment_id = experiment.experiment_id

    logger.info(
        f"Using MLflow experiment ID: {experiment_id}"
    )

    results = []

    for model_index, model_config in enumerate(
        models,
        start=1,
    ):
        model_id = model_config.get(
            "model_id",
            f"unnamed_model_{model_index}",
        )

        logger.info("=" * 80)
        logger.info(
            f"[{model_index}/{len(models)}] "
            f"Migrating model: {model_id}"
        )
        logger.info("=" * 80)

        try:
            result = _migrate_single_model(
                model_config=model_config,
                experiment_id=experiment_id,
                metrics_config=metrics_config,
                skip_existing=skip_existing,
                calculate_checksums=calculate_checksums,
                require_optional_artifacts=require_optional_artifacts,
            )

        except Exception as error:
            logger.exception(
                f"Failed to migrate model '{model_id}'"
            )

            result = {
                "model_id": model_id,
                "status": "failed",
                "run_id": None,
                "checkpoint_path": model_config.get(
                    "checkpoint_path"
                ),
                "metrics_csv_path": model_config.get(
                    "metrics_csv_path"
                ),
                "checkpoint_sha256": None,
                "n_tasks": None,
                "n_per_task_metrics": None,
                "n_summary_metrics": None,
                "error": (
                    f"{type(error).__name__}: {error}"
                ),
            }

            results.append(result)

            if not continue_on_error:
                logger.error(
                    "Stopping migration because "
                    "continue_on_error=false"
                )
                break

            continue

        results.append(result)

    results_df = pd.DataFrame(results)

    report_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    results_df.to_csv(
        report_path,
        index=False,
    )

    status_counts = (
        results_df["status"]
        .value_counts(dropna=False)
        .to_dict()
    )

    n_completed = int(
        status_counts.get("completed", 0)
    )

    n_skipped = int(
        status_counts.get("skipped_existing", 0)
    )

    n_failed = int(
        status_counts.get("failed", 0)
    )

    logger.info("=" * 80)
    logger.info("Historical ADME model migration finished")
    logger.info(f"Completed: {n_completed}")
    logger.info(f"Skipped existing: {n_skipped}")
    logger.info(f"Failed: {n_failed}")
    logger.info(f"Migration report: {report_path}")
    logger.info("=" * 80)

    context["mlflow_migration_results"] = results_df
    context["mlflow_migration_report_path"] = str(
        report_path
    )
    context["mlflow_experiment_id"] = experiment_id
    context["mlflow_tracking_uri"] = tracking_uri

    return {
        "mlflow_migration_results": results_df,
        "mlflow_migration_report_path": str(report_path),
        "mlflow_experiment_id": experiment_id,
        "mlflow_tracking_uri": tracking_uri,
    }