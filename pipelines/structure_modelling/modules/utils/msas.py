import subprocess
import shlex
import hashlib
import shutil

from pathlib import Path

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

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

    logger.info(
        "[MSA] Running: " + " ".join(cmd)
    )

    subprocess.run(
        cmd,
        check=True,
    )

    generated_a3m = find_generated_a3m(
        search_out_dir,
        query_stem=query_stem,
    )

    if generated_a3m is None:
        raise FileNotFoundError(
            f"No A3M file was generated for {sequence_id} under {search_out_dir}"
        )

    shutil.copy(
        generated_a3m,
        final_a3m,
    )

    logger.info(
        f"[MSA] Cached MSA for {sequence_id}: {final_a3m}"
    )

    return final_a3m
