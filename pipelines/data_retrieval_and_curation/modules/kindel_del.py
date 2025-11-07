# modules/kindel_data.py
from pipeline.task_registry import register_task
import s3fs
from tqdm import tqdm
import threading
import os
import time
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_task(
    "retrieve_kindel_datasets",
    category="DNA-Encoded Libraries",
    description="Download multiple KinDel datasets from public S3."
)
def retrieve_kindel_datasets(config, data=None):
    """
    Downloads multiple datasets from S3 with resumable support.

    Config keys:
        datasets: list of dicts with 'name' and 'file'
        s3_prefix: common S3 folder/prefix
        output: dict with 'directory' and 'overwrite'
    """
    datasets = config.get("datasets", [])
    s3_prefix = config.get("s3_prefix")
    output_cfg = config.get("output", {})
    directory = output_cfg.get("directory", ".")
    overwrite = output_cfg.get("overwrite", False)

    if not datasets:
        raise ValueError("Config must specify 'datasets' list.")
    if not s3_prefix:
        raise ValueError("Config must specify 's3_prefix'.")

    os.makedirs(directory, exist_ok=True)
    fs = s3fs.S3FileSystem(anon=True)
    results = {}

    for ds in datasets:
        name = ds["name"]
        filename = ds["file"]
        remote = f"{s3_prefix}/{filename}"
        local = os.path.join(directory, filename)

        # Check file existence
        size = fs.info(remote)["size"]
        existing_size = 0
        if os.path.exists(local):
            existing_size = os.path.getsize(local)
            if existing_size >= size:
                if overwrite:
                    logger.info(f"Overwriting existing file: {local}")
                    existing_size = 0
                else:
                    logger.info(f"File already exists: {local}")
                    results[name] = local
                    continue

        logger.info(f"Starting download: {name} ({size/1e6:.2f} MB) -> {local}")
        if existing_size > 0:
            logger.info(f"Resuming from {existing_size/1e6:.2f} MB")

        # Progress bar
        pbar = tqdm(total=size, unit="B", unit_scale=True, desc=f"{name} Download", initial=existing_size)

        # Monitor thread
        def monitor_file(path, total_size, bar):
            while True:
                if os.path.exists(path):
                    current_size = os.path.getsize(path)
                    bar.update(current_size - bar.n)
                    if current_size >= total_size:
                        break
                time.sleep(0.2)
            bar.close()

        t = threading.Thread(target=monitor_file, args=(local, size, pbar))
        t.start()

        # Download with resume
        with fs.open(remote, "rb") as s3f, open(local, "ab") as f_local:
            if existing_size > 0:
                s3f.seek(existing_size)
            chunk_size = 1024 * 1024
            while True:
                chunk = s3f.read(chunk_size)
                if not chunk:
                    break
                f_local.write(chunk)

        t.join()
        logger.info(f"Download complete: {local}")
        results[name] = local

    return {"datasets": results}
