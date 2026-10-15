import argparse
import os
import sys
import traceback
import yaml

from pathlib import Path

from pipeline.logger import setup_logger


logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)


SCHRODINGER_BACKENDS = {
    "schrodingerglideem",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run docking pipeline via dynamically selected workflow."
        )
    )

    parser.add_argument(
        "-p",
        "--params",
        required=True,
        help="Path to the YAML configuration file."
    )

    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug output."
    )

    return parser.parse_args()


def load_config(config_path: Path):

    if not config_path.is_file():
        logger.error(
            f"Config file not found: {config_path}"
        )
        sys.exit(1)

    with config_path.open("r") as handle:
        config = yaml.safe_load(handle)

    if not isinstance(config, dict):
        logger.error(
            f"Invalid YAML configuration: {config_path}"
        )
        sys.exit(1)

    return config


def requires_schrodinger(config: dict) -> bool:

    backend_name = (
        config
        .get("backend", {})
        .get("name", "")
        .lower()
    )

    return backend_name in SCHRODINGER_BACKENDS


def running_under_schrodinger() -> bool:

    return (
        os.environ.get("KITBAG_SCHRODINGER_PYTHON")
        == "1"
    )


def relaunch_with_schrodinger():

    schrodinger = os.environ.get("SCHRODINGER")

    if not schrodinger:
        logger.error(
            "This workflow requires Schrodinger, "
            "but $SCHRODINGER is not set."
        )
        sys.exit(1)

    schrodinger_run = Path(schrodinger) / "run"

    if not schrodinger_run.is_file():
        logger.error(
            f"Schrodinger launcher not found: "
            f"{schrodinger_run}"
        )
        sys.exit(1)

    script_path = Path(__file__).resolve()

    cmd = [
        str(schrodinger_run),
        "python3",
        str(script_path),
        *sys.argv[1:]
    ]

    logger.info(
        "Schrodinger backend selected. "
        "Relaunching workflow under Schrodinger Python."
    )

    if os.environ.get("KITBAG_DEBUG_LAUNCH") == "1":
        logger.info(
            f"Launch command: {' '.join(cmd)}"
        )

    env = os.environ.copy()

    env["KITBAG_SCHRODINGER_PYTHON"] = "1"

    os.execve(
        str(schrodinger_run),
        cmd,
        env
    )


def main():

    args = parse_args()

    config_path = (
        Path(args.params)
        .expanduser()
        .resolve()
    )

    # Load the config before importing workflows.
    config = load_config(config_path)

    # --------------------------------------------------
    # Select Python runtime
    # --------------------------------------------------

    if (
        requires_schrodinger(config)
        and not running_under_schrodinger()
    ):
        relaunch_with_schrodinger()

    # --------------------------------------------------
    # Workflow discovery happens only after the correct
    # Python runtime has been selected.
    # --------------------------------------------------

    from workflows import (
        get_workflow,
        list_workflows,
    )

    workflow_name = config.get("workflow_name")

    if not workflow_name:
        logger.error(
            "Config must include a 'workflow_name' key."
        )
        sys.exit(1)

    workflow_func = get_workflow(workflow_name)

    if workflow_func is None:

        logger.error(
            f"Unknown workflow: '{workflow_name}'"
        )

        logger.info(
            "Available workflows:"
        )

        for name in list_workflows():
            logger.info(
                f"  - {name}"
            )

        sys.exit(1)

    if args.debug:

        backend_name = (
            config
            .get("backend", {})
            .get("name")
        )

        logger.info(
            f"Selected workflow: {workflow_name}"
        )

        logger.info(
            f"Backend: {backend_name}"
        )

        logger.info(
            f"Config path: {config_path}"
        )

        runtime = (
            "Schrodinger Python"
            if running_under_schrodinger()
            else "standard Python"
        )

        logger.info(
            f"Runtime: {runtime}"
        )

    # --------------------------------------------------
    # Run workflow
    # --------------------------------------------------

    try:

        workflow_func(
            str(config_path)
        )

    except KeyboardInterrupt:

        logger.warning(
            f"Workflow '{workflow_name}' "
            "interrupted by user."
        )

        sys.exit(130)

    except Exception as exc:

        logger.error(
            f"[FATAL] Workflow "
            f"'{workflow_name}' failed: {exc}"
        )

        traceback.print_exc()

        sys.exit(1)


if __name__ == "__main__":
    main()