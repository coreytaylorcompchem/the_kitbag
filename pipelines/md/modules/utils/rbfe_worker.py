import os
import yaml
import time
import subprocess
import glob

from pathlib import Path
import multiprocessing as mp
from queue import Empty

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def _get_openfe_progress(workdir):

    matches = glob.glob(
        os.path.join(
            workdir,
            "**",
            "simulation_real_time_analysis.yaml"
        ),
        recursive=True
    )

    if not matches:
        return None

    try:

        with open(matches[0]) as f:
            data = yaml.safe_load(f)

        if not data:
            return None

        latest = data[-1]

        return latest.get(
            "percent_complete"
        )

    except Exception:

        return None

def _run_single_transformation(
    transform_json,
    output_dir,
    env,
):

    stem = Path(transform_json).stem

    result_json = os.path.join(
        output_dir,
        f"{stem}.json",
    )

    working_dir = os.path.join(
        output_dir,
        stem,
    )

    #
    # Skip completed calculations
    #

    if os.path.exists(result_json):

        logger.info(
            f"Skipping completed {stem}"
        )

        return True

    logger.info(
        f"Launching on visible GPU {env['CUDA_VISIBLE_DEVICES']}"
    )
    
    logger.info(
        f"Running {stem}"
    )

    cmd = [
        "bash",
        "-c",
        (
            "openfe quickrun "
            f"{transform_json} "
            f"-o {result_json} "
            f"-d {working_dir}"
        )
    ]

    proc = subprocess.Popen(
        cmd,
        env=env,
    )

    last_report = -1

    while proc.poll() is None:

        progress = _get_openfe_progress(
            working_dir
        )

        if (
            progress is not None
            and int(progress) != last_report
        ):

            logger.info(
                f"[GPU {env['CUDA_VISIBLE_DEVICES']}] "
                f"{stem}: "
                f"{progress:.1f}%"
            )

            last_report = int(progress)

        time.sleep(300)

    if proc.returncode != 0:

        logger.error(
            f"FAILED {stem}"
        )

        return False

    logger.info(
        f"Completed {stem}"
    )

    return True

def _openfe_worker(
    gpu_id,
    queue,
    failed,
    output_dir,
):

    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)

    while not queue.empty():

        try:
            transform_json = queue.get_nowait()

        except Exception:
            break

        ok = _run_single_transformation(
            transform_json,
            output_dir,
            env,
        )

        if not ok:

            failed.append(
                Path(transform_json).stem
            )