import os
import psutil
import time
import threading

import subprocess
import shlex
import hashlib
import shutil

from pathlib import Path

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def directory_summary(path: Path):
    if not path.exists():
        return {
            "files": 0,
            "size_mb": 0.0,
        }

    total_size = 0
    file_count = 0

    for f in path.rglob("*"):
        if f.is_file():
            file_count += 1
            try:
                total_size += f.stat().st_size
            except OSError:
                pass

    return {
        "files": file_count,
        "size_mb": total_size / (1024 * 1024),
    }


def run_subprocess_streaming(
    cmd,
    cwd: Path = None,
    monitor_dir: Path = None,
    log_prefix: str = "[CMD]",
    heartbeat_seconds: int = 10,
):
    """
    Run an external command while streaming stdout/stderr into the pipeline logger.

    Also emits periodic heartbeat messages so long-running MSA searches do not
    look stalled.
    """

    start_time = time.time()
    stop_event = threading.Event()

    def heartbeat():
        while not stop_event.wait(heartbeat_seconds):
            elapsed = time.time() - start_time

            if monitor_dir is not None:
                summary = directory_summary(monitor_dir)

                mem = psutil.virtual_memory()

                logger.info(
                    f"{log_prefix} RUNTIME MSA: {elapsed / 60:.1f} min "
                    f"| files={summary['files']} "
                    f"| size={summary['size_mb']:.1f} MB "
                    # f"| monitor_dir={monitor_dir}"
                )

                logger.debug(
                    f"{log_prefix} RAM: "
                    f"used={mem.used/1024**3:.1f}GB "
                    f"avail={mem.available/1024**3:.1f}GB"
                )
            else:
                logger.info(
                    f"{log_prefix} still running after {elapsed / 60:.1f} min"
                )

    heartbeat_thread = threading.Thread(
        target=heartbeat,
        daemon=True,
    )

    heartbeat_thread.start()

    logger.info(
        f"{log_prefix} === MSA STEP STARTED ==="
    )

    process = subprocess.Popen(
        cmd,
        cwd=str(cwd) if cwd is not None else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        universal_newlines=True,
    )

    try:
        for line in process.stdout:
            line = line.rstrip()
            if line:
                logger.info(f"{log_prefix} {line}")

        return_code = process.wait()

    finally:
        stop_event.set()
        heartbeat_thread.join(timeout=5)

    elapsed = time.time() - start_time

    if return_code != 0:
        raise subprocess.CalledProcessError(
            return_code,
            cmd,
        )

    logger.info(
        f"{log_prefix} command completed successfully in {elapsed / 60:.1f} min"
    )

def sequence_hash(sequence: str) -> str:
    sequence = sequence.strip().upper()
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()[:20]


def safe_name(value: str) -> str:
    return "".join(
        c if c.isalnum() or c in ("-", "_") else "_"
        for c in str(value)
    )


def command_to_list(command):
    if isinstance(command, list):
        return [str(x) for x in command]

    if isinstance(command, str):
        return shlex.split(command)

    raise ValueError(
        "msa.command must be either a string or a list"
    )


def write_fasta(sequence_id: str, sequence: str, fasta_path: Path):
    fasta_path.parent.mkdir(parents=True, exist_ok=True)

    sequence = sequence.strip().upper()

    with open(fasta_path, "w") as f:
        f.write(f">{sequence_id}\n")

        for i in range(0, len(sequence), 80):
            f.write(sequence[i:i + 80] + "\n")


def find_generated_a3m(search_dir: Path, query_stem: str):
    candidates = sorted(search_dir.rglob("*.a3m"))

    if not candidates:
        return None

    exact = [
        p for p in candidates
        if p.stem == query_stem
    ]

    if exact:
        return exact[0]

    containing = [
        p for p in candidates
        if query_stem in p.stem
    ]

    if containing:
        return containing[0]

    candidates = sorted(
        candidates,
        key=lambda p: p.stat().st_size,
        reverse=True
    )

    return candidates[0]


def generate_local_msa_for_sequence(
    sequence_id: str,
    sequence: str,
    database_dir: Path,
    cache_dir: Path,
    command,
    extra_args=None,
    reuse_existing=True,
    overwrite=False,
):
    seq_hash = sequence_hash(sequence)

    a3m_cache_dir = cache_dir / "a3m"
    work_root = cache_dir / "work"

    a3m_cache_dir.mkdir(parents=True, exist_ok=True)
    work_root.mkdir(parents=True, exist_ok=True)

    final_a3m = a3m_cache_dir / f"{seq_hash}.a3m"

    if final_a3m.exists() and reuse_existing and not overwrite:
        logger.info(
            f"[MSA] Reusing cached MSA for {sequence_id}: {final_a3m}"
        )
        return final_a3m

    work_dir = work_root / seq_hash

    if work_dir.exists() and overwrite:
        shutil.rmtree(work_dir)

    work_dir.mkdir(parents=True, exist_ok=True)

    query_stem = safe_name(sequence_id)
    fasta_path = work_dir / f"{query_stem}.fasta"
    search_out_dir = work_dir / "search"

    if search_out_dir.exists() and overwrite:
        shutil.rmtree(search_out_dir)

    search_out_dir.mkdir(parents=True, exist_ok=True)

    write_fasta(
        sequence_id=sequence_id,
        sequence=sequence,
        fasta_path=fasta_path,
    )

    cmd = (
        command_to_list(command)
        + [
            str(fasta_path),
            str(database_dir),
            str(search_out_dir),
        ]
        + [str(x) for x in (extra_args or [])]
    )

    if not database_dir.exists():
        raise FileNotFoundError(
            f"MSA database directory does not exist: {database_dir}"
        )

    if shutil.which(cmd[0]) is None:
        raise FileNotFoundError(
            f"MSA command not found: {cmd[0]}\n"
            f"Full command was: {' '.join(cmd)}\n"
            "If using a separate MSA environment, set for example:\n"
            "  command: \"conda run -n msa_tools colabfold_search\""
        )
    
    logger.debug(
    f"[MSA] command executable: {cmd[0]}"
    )

    logger.debug(
        f"[MSA] database_dir: {database_dir}"
    )

    logger.debug(
        f"[MSA] work_dir: {work_dir}"
    )

    logger.debug(
        f"[MSA] search_out_dir: {search_out_dir}"
    )

    logger.debug(
        f"[MSA] final_a3m: {final_a3m}"
    )

    logger.info(
        "[MSA] Running: " + " ".join(cmd)
    )

    run_subprocess_streaming(
        cmd=cmd,
        cwd=work_dir,
        monitor_dir=search_out_dir,
        log_prefix=f"[MSA:{sequence_id}]",
    )

    generated_a3m = find_generated_a3m(
        search_out_dir,
        query_stem=query_stem,
    )

    if generated_a3m is None:
        produced_files = sorted(
            str(p.relative_to(search_out_dir))
            for p in search_out_dir.rglob("*")
            if p.is_file()
        )

        preview = produced_files[:50]

        raise FileNotFoundError(
            f"No A3M file was generated for {sequence_id} under {search_out_dir}\n"
            f"Number of produced files: {len(produced_files)}\n"
            f"First files:\n"
            + "\n".join(preview)
        )

    shutil.copy(
        generated_a3m,
        final_a3m,
    )

    logger.info(
        f"[MSA] Cached MSA for {sequence_id}: {final_a3m}"
    )

    return final_a3m
