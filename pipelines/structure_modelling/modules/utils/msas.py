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

def detect_colabfold_stage(search_dir: Path):
    """
    Infer progress through the ColabFold/MMseqs2 monomer MSA workflow
    from intermediate database files.

    Stage numbers are operational milestones only. They must not be
    interpreted as percentage completion because stages can differ
    substantially in runtime.
    """

    search_dir = Path(search_dir)

    if not search_dir.exists():
        return 1, 10, "PREPARING_QUERIES"

    def exists(prefix):
        """
        Check whether an MMseqs2 database/output with this prefix exists.

        MMseqs2 databases consist of multiple associated files, commonly
        including a .dbtype file, so check the prefix and related files.
        """
        return (
            (search_dir / prefix).exists()
            or (search_dir / f"{prefix}.dbtype").exists()
            or any(search_dir.glob(f"{prefix}.*"))
        )

    # --------------------------------------------------------------
    # 10. Final per-query A3Ms are being written
    # --------------------------------------------------------------

    if any(search_dir.glob("*.a3m")):
        return 10, 10, "WRITING_A3MS"

    # --------------------------------------------------------------
    # 9. UniRef and environmental MSAs have been combined
    # --------------------------------------------------------------

    if exists("final.a3m"):
        return 9, 10, "MERGING_MSAS"

    # --------------------------------------------------------------
    # 8. Environmental alignments are being filtered / converted
    #    into the environmental MSA
    # --------------------------------------------------------------

    if (
        exists("bfd.mgnify30.metaeuk30.smag30.a3m")
        or exists("res_env_exp_realign_filter")
        or exists("res_env_exp_realign")
    ):
        return 8, 10, "BUILDING_ENV_MSA"

    # --------------------------------------------------------------
    # 7. Environmental database hits are being expanded
    # --------------------------------------------------------------

    if exists("res_env_exp"):
        return 7, 10, "EXPANDING_ENV"

    # --------------------------------------------------------------
    # 6. Searching environmental sequence database
    # --------------------------------------------------------------

    if exists("res_env"):
        return 6, 10, "SEARCHING_ENV"

    # --------------------------------------------------------------
    # 5. UniRef alignments are being filtered / converted
    #    into the UniRef MSA
    # --------------------------------------------------------------

    if (
        exists("uniref.a3m")
        or exists("res_exp_realign_filter")
    ):
        return 5, 10, "BUILDING_UNIREF_MSA"

    # --------------------------------------------------------------
    # 4. Detailed realignment of expanded UniRef hits
    # --------------------------------------------------------------

    if exists("res_exp_realign"):
        return 4, 10, "ALIGNING_UNIREF"

    # --------------------------------------------------------------
    # 3. Initial UniRef hits are being expanded
    # --------------------------------------------------------------

    if exists("res_exp"):
        return 3, 10, "EXPANDING_UNIREF"

    # --------------------------------------------------------------
    # 2. Searching UniRef
    # --------------------------------------------------------------

    if exists("res"):
        return 2, 10, "SEARCHING_UNIREF"

    # --------------------------------------------------------------
    # 1. FASTA has been converted / is being converted into qdb
    # --------------------------------------------------------------

    if exists("qdb"):
        return 1, 10, "PREPARING_QUERIES"

    return 1, 10, "PREPARING_QUERIES"

def run_subprocess_streaming(
    cmd,
    cwd: Path = None,
    monitor_dir: Path = None,
    log_prefix: str = "[CMD]",
    heartbeat_seconds: int = 60,
    expected_queries: int = None,
):
    """
    Run an external command while streaming stdout/stderr into the pipeline logger.

    Also emits periodic heartbeat messages so long-running MSA searches do not
    look stalled.

    For ColabFold MSA searches, the heartbeat reports:
      - total runtime
      - inferred MMseqs2/ColabFold stage
      - time spent in the current stage
      - number of final A3M files written
      - files and disk usage in the search directory

    Stage numbers are operational milestones only. They do not represent equal
    fractions of total runtime.
    """

    start_time = time.time()
    stop_event = threading.Event()

    progress_state = {
        "stage_key": None,
        "stage_started": start_time,
    }

    def heartbeat():
        while not stop_event.wait(heartbeat_seconds):
            elapsed = time.time() - start_time

            if monitor_dir is not None:
                summary = directory_summary(monitor_dir)

                stage_num, stage_total, stage_name = detect_colabfold_stage(
                    monitor_dir
                )

                stage_key = (stage_num, stage_name)

                if progress_state["stage_key"] != stage_key:
                    progress_state["stage_key"] = stage_key
                    progress_state["stage_started"] = time.time()

                    logger.info(
                        f"{log_prefix} Stage changed to "
                        f"{stage_num}/{stage_total} ({stage_name})"
                    )

                stage_elapsed = (
                    time.time() - progress_state["stage_started"]
                )

                a3m_count = sum(
                    1
                    for p in monitor_dir.rglob("*.a3m")
                    if p.is_file()
                )

                mem = psutil.virtual_memory()

                if expected_queries is not None:
                    a3m_text = f"{a3m_count}/{expected_queries}"
                else:
                    a3m_text = str(a3m_count)

                logger.info(
                    f"{log_prefix} RUNTIME MSA: "
                    f"{elapsed / 60:.1f} min "
                    f"| stage={stage_num}/{stage_total} ({stage_name}) "
                    f"| stage time elapsed={stage_elapsed / 60:.1f} min "
                    f"| A3Ms DONE={a3m_text} "
                    f"| files={summary['files']} "
                    # f"| size={summary['size_mb'\]:.1f} MB"
                )

                logger.debug(
                    f"{log_prefix} RAM: "
                    f"used={mem.used/1024**3:.1f}GB "
                    f"avail={mem.available/1024**3:.1f}GB"
                )

            else:
                logger.info(
                    f"{log_prefix} still running after "
                    f"{elapsed / 60:.1f} min"
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
        f"{log_prefix} command completed successfully in "
        f"{elapsed / 60:.1f} min"
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

def write_multi_fasta(sequence_jobs, fasta_path: Path):
    """
    Write multiple unique sequences to a single FASTA file.

    Parameters
    ----------
    sequence_jobs
        Iterable of dictionaries containing:
            sequence_id
            sequence
    fasta_path
        Output FASTA path.

    Notes
    -----
    sequence_id should be filesystem-safe because ColabFold commonly uses
    the FASTA identifier when naming the resulting A3M file.
    """
    fasta_path.parent.mkdir(parents=True, exist_ok=True)

    seen_ids = set()

    with open(fasta_path, "w") as f:
        for job in sequence_jobs:
            sequence_id = safe_name(job["sequence_id"])
            sequence = job["sequence"].strip().upper()

            if not sequence_id:
                raise ValueError("Empty sequence_id in batch MSA input")

            if sequence_id in seen_ids:
                raise ValueError(
                    f"Duplicate sequence_id in batch MSA input: {sequence_id}"
                )

            if not sequence:
                raise ValueError(
                    f"Empty protein sequence for sequence_id={sequence_id}"
                )

            seen_ids.add(sequence_id)

            f.write(f">{sequence_id}\n")

            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i + 80] + "\n")

def generate_local_msas_batch(
    sequence_jobs,
    database_dir: Path,
    cache_dir: Path,
    command,
    extra_args=None,
    reuse_existing=True,
    overwrite=False,
):
    """
    Generate MSAs for multiple unique sequences in one colabfold_search call.

    Parameters
    ----------
    sequence_jobs
        Iterable of dictionaries containing:
            sequence_id: stable filesystem-safe identifier
            sequence: protein sequence

        In this pipeline, sequence_id should be sequence_hash(sequence).

    database_dir
        Local ColabFold/MMseqs2 database directory.

    cache_dir
        Root MSA cache directory.

    command
        Command string or list, for example:
            "colabfold_search"
        or:
            "conda run -n msa_tools colabfold_search"

    extra_args
        Additional arguments passed to colabfold_search.

    reuse_existing
        Reuse existing <cache_dir>/a3m/<sequence_hash>.a3m files.

    overwrite
        Regenerate MSAs even when cache files already exist.

    Returns
    -------
    dict
        Mapping from sequence_id to cached A3M Path.
    """
    sequence_jobs = list(sequence_jobs)

    if not sequence_jobs:
        logger.info("[MSA] No sequence jobs supplied to batch generator")
        return {}

    database_dir = Path(database_dir)
    cache_dir = Path(cache_dir)

    a3m_cache_dir = cache_dir / "a3m"
    work_root = cache_dir / "work"
    batch_work_dir = work_root / "batch"
    search_out_dir = batch_work_dir / "search"
    fasta_path = batch_work_dir / "queries.fasta"

    a3m_cache_dir.mkdir(parents=True, exist_ok=True)
    work_root.mkdir(parents=True, exist_ok=True)

    if not database_dir.exists():
        raise FileNotFoundError(
            f"MSA database directory does not exist: {database_dir}"
        )

    # --------------------------------------------------------------
    # Determine which MSAs can be reused and which require searching.
    # --------------------------------------------------------------

    msa_paths = {}
    pending_jobs = []

    seen_sequence_ids = set()

    for job in sequence_jobs:
        sequence_id = safe_name(job["sequence_id"])
        sequence = job["sequence"].strip().upper()

        if sequence_id in seen_sequence_ids:
            raise ValueError(
                f"Duplicate batch MSA sequence_id: {sequence_id}"
            )

        seen_sequence_ids.add(sequence_id)

        expected_hash = sequence_hash(sequence)

        if sequence_id != expected_hash:
            raise ValueError(
                f"Batch MSA sequence_id does not match sequence hash: "
                f"sequence_id={sequence_id}, expected={expected_hash}"
            )

        final_a3m = a3m_cache_dir / f"{sequence_id}.a3m"

        if (
            final_a3m.exists()
            and reuse_existing
            and not overwrite
        ):
            logger.info(
                f"[MSA] Reusing cached MSA for {sequence_id}: "
                f"{final_a3m}"
            )
            msa_paths[sequence_id] = final_a3m
        else:
            pending_jobs.append({
                "sequence_id": sequence_id,
                "sequence": sequence,
            })

    logger.info(
        f"[MSA] Batch cache lookup complete | "
        f"total={len(sequence_jobs)} "
        f"| cached={len(msa_paths)} "
        f"| pending={len(pending_jobs)}"
    )

    if not pending_jobs:
        logger.info("[MSA] All requested MSAs were found in cache")
        return msa_paths

    # --------------------------------------------------------------
    # Use a clean batch work directory for the pending sequences.
    #
    # This prevents an old A3M from being mistaken for output from the
    # current search.
    # --------------------------------------------------------------

    if batch_work_dir.exists():
        shutil.rmtree(batch_work_dir)

    search_out_dir.mkdir(parents=True, exist_ok=True)

    write_multi_fasta(
        sequence_jobs=pending_jobs,
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

    if shutil.which(cmd[0]) is None:
        raise FileNotFoundError(
            f"MSA command not found: {cmd[0]}\n"
            f"Full command was: {' '.join(cmd)}\n"
            "If using a separate MSA environment, set for example:\n"
            '  command: "conda run -n msa_tools colabfold_search"'
        )

    logger.debug(f"[MSA] command executable: {cmd[0]}")
    logger.debug(f"[MSA] database_dir: {database_dir}")
    logger.debug(f"[MSA] batch_work_dir: {batch_work_dir}")
    logger.debug(f"[MSA] search_out_dir: {search_out_dir}")
    logger.debug(f"[MSA] batch_fasta: {fasta_path}")

    logger.info(
        f"[MSA] Running one batched search for "
        f"{len(pending_jobs)} unique sequences"
    )

    logger.info(
        "[MSA] Running: " + " ".join(cmd)
    )

    run_subprocess_streaming(
        cmd=cmd,
        cwd=batch_work_dir,
        monitor_dir=search_out_dir,
        log_prefix="[MSA:BATCH]",
        expected_queries=len(pending_jobs),
    )

    # --------------------------------------------------------------
    # Discover and validate the generated A3M files.
    # --------------------------------------------------------------

    generated_a3ms = sorted(search_out_dir.rglob("*.a3m"))

    generated_by_stem = {}

    for path in generated_a3ms:
        if path.stem in generated_by_stem:
            raise RuntimeError(
                f"Multiple A3M files generated with stem '{path.stem}': "
                f"{generated_by_stem[path.stem]} and {path}"
            )

        generated_by_stem[path.stem] = path

    missing_sequence_ids = []

    for job in pending_jobs:
        sequence_id = job["sequence_id"]
        generated_a3m = generated_by_stem.get(sequence_id)

        if generated_a3m is None:
            missing_sequence_ids.append(sequence_id)
            continue

        final_a3m = a3m_cache_dir / f"{sequence_id}.a3m"

        shutil.copy2(
            generated_a3m,
            final_a3m,
        )

        msa_paths[sequence_id] = final_a3m

        logger.info(
            f"[MSA] Cached batch MSA for {sequence_id}: "
            f"{final_a3m}"
        )

    if missing_sequence_ids:
        produced_files = [
            str(path.relative_to(search_out_dir))
            for path in generated_a3ms
        ]

        raise FileNotFoundError(
            "The batched MSA search completed, but some expected A3M "
            "files were not produced.\n"
            f"Missing sequence IDs ({len(missing_sequence_ids)}):\n"
            + "\n".join(missing_sequence_ids[:50])
            + "\n\nGenerated A3M files:\n"
            + "\n".join(produced_files[:100])
        )

    if len(msa_paths) != len(sequence_jobs):
        raise RuntimeError(
            f"MSA batch result count mismatch: "
            f"requested={len(sequence_jobs)}, resolved={len(msa_paths)}"
        )

    logger.info(
        f"[MSA] Batch MSA generation complete | "
        f"resolved={len(msa_paths)}"
    )

    return msa_paths