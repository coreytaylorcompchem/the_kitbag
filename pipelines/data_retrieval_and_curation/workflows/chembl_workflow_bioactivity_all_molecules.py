from workflows import register_workflow
from chembl_webresource_client.new_client import new_client

from pipeline.parallel_runner import ParallelWorkflowRunner
from pipeline.task_registry import get_task
from pipeline.logger import setup_logger

import math
import os
import json
import requests
import time
from math import ceil
import pandas as pd

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

import requests
import time
from tqdm import tqdm

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def fetch_all_molecule_ids_webclient(config):

    molecule_client = new_client.molecule
    max_records = config.get("max_records", 100000)
    if isinstance(max_records, str) and max_records.lower() == "all":
        max_records = None  # Fetch all!

    resume = config.get("resume", False)
    cleanup = config.get("cleanup", False)
    output_dir = config.get("output", {}).get("directory", "outputs/chembl_molecules")
    checkpoint_file = os.path.join(output_dir, "checkpoint.json")

    os.makedirs(output_dir, exist_ok=True)

    all_ids = []
    offset = 0

    # Set up resume, reading in json file if previous queries were interrupted.
    if resume and os.path.exists(checkpoint_file):
        with open(checkpoint_file, "r") as f:
            checkpoint = json.load(f)
            all_ids = checkpoint.get("collected_ids", [])
            offset = checkpoint.get("last_index", 0)
        logger.info(f"Resuming from offset {offset}, {len(all_ids)} IDs loaded")

    limit = 1000  # max page size supported by ChEMBL
    pbar = tqdm(desc="Fetching molecules", initial=len(all_ids), total=max_records)

    def fetch_page(offset, limit):
        return molecule_client.filter().only("molecule_chembl_id")[offset:offset + limit]

    while True:
        try:
            page = retry_fetch(lambda: fetch_page(offset, limit))
            if not page:
                logger.info("No more results. Finished.")
                break

            for mol in page:
                if "molecule_chembl_id" in mol:
                    all_ids.append(mol["molecule_chembl_id"])
                    offset += 1 
                    pbar.update(1)

                    if resume and offset % 100 == 0:
                        with open(checkpoint_file, "w") as f:
                            json.dump({"last_index": offset, "collected_ids": all_ids}, f)

                if max_records and len(all_ids) >= max_records:
                    logger.info("Reached max_records limit.")
                    raise StopIteration

            if resume and offset % 1000 == 0:
                with open(checkpoint_file, "w") as f:
                    json.dump({"last_index": offset, "collected_ids": all_ids}, f)

        except StopIteration:
            break
        except Exception as e:
            logger.warning(f"Error while fetching: {e}. Retrying...")
            time.sleep(10)

    pbar.close()

    if resume:
        with open(checkpoint_file, "w") as f:
            json.dump({"last_index": offset, "collected_ids": all_ids}, f)

    if resume and cleanup and os.path.exists(checkpoint_file):
        try:
            os.remove(checkpoint_file)
            logger.info("Checkpoint file deleted.")
        except Exception as e:
            logger.warning(f"Failed to delete checkpoint file: {e}")

    logger.info(f"Total molecule IDs collected: {len(all_ids)}")
    return all_ids


def retry_fetch(fetch_fn, retries=5, backoff=5):
    import random
    for i in range(retries):
        try:
            return fetch_fn()
        except requests.exceptions.RequestException as e:
            wait = backoff * (2 ** i) + random.uniform(0, 1)
            logger.warning(f"Retry {i + 1}/{retries} after error: {e}. Waiting {wait:.1f}s...")
            time.sleep(wait)
    raise RuntimeError("❌ Maximum retries exceeded.")


@register_workflow("chembl_all_molecules",
                   description="Download ChEMBL molecule metadata in parallel.")
def run_chembl_all_molecules_parallel(config):
    max_records = config.get("max_records", 100000)
    batch_size = config.get("batch_size", 1000)

    logger.info(f"Retrieving up to {max_records} ChEMBL molecule IDs...")
    
    all_ids = fetch_all_molecule_ids_webclient(config)

    logger.info(f"Retrieved {len(all_ids)} molecule IDs. Preparing batches...")

    def chunk_list(lst, chunk_size):
        for i in range(0, len(lst), chunk_size):
            yield lst[i:i + chunk_size]

    batches = list(chunk_list(all_ids, batch_size))

    runner = ParallelWorkflowRunner(
        workflow_func=_single_batch_runner,
        config={**config, "molecule_batches": [(i, batch) for i, batch in enumerate(batches)]},
        input_key="molecule_batches",
        output_key="batch_index",
        filename_pattern="batch_{batch_index}.csv",
        combined_filename=config.get("output", {}).get("filename", "all_molecules_combined.csv"),
        output_dir=config.get("output", {}).get("directory", "outputs/chembl_parallel"),
        input_is_pair=True,
        use_multiprocessing=True
    )

    runner.inputs = [(i, batch) for i, batch in enumerate(batches)]
    return runner.run()

def _single_batch_runner(input_pair):
    batch_index, chembl_ids = input_pair
    task = get_task("retrieve_compound_metadata_batch")
    if not task:
        raise RuntimeError("Task 'retrieve_compound_metadata_batch' not registered.")
    local_config = {"molecule_chembl_ids": chembl_ids}
    return task(local_config)