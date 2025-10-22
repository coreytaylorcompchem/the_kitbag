from workflows import register_workflow
from chembl_webresource_client.new_client import new_client

from pipeline.parallel_runner import ParallelWorkflowRunner
from pipeline.task_registry import get_task
from pipeline.logger import setup_logger

import math
import random
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
        max_records = None  # Fetch all

    resume = config.get("resume", False)
    cleanup = config.get("cleanup", False)
    output_dir = config.get("output", {}).get("directory", "outputs/chembl_molecules")
    checkpoint_file = os.path.join(output_dir, "checkpoint.json")

    os.makedirs(output_dir, exist_ok=True)

    all_ids = []
    current_index = 0

    if resume and os.path.exists(checkpoint_file):
        with open(checkpoint_file, "r") as f:
            checkpoint = json.load(f)
            all_ids = checkpoint.get("collected_ids", [])
            current_index = checkpoint.get("last_index", 0)
        logger.info(f"⏸ Resuming from index {current_index}, {len(all_ids)} IDs loaded")

    logger.info(f"Fetching up to {max_records if max_records else 'all'} molecule IDs...")

    max_retries = 5
    retry_wait = 1  # base seconds

    pbar = tqdm(total=max_records, desc="Fetching molecules", initial=len(all_ids))

    while not max_records or len(all_ids) < max_records:
        attempt = 0
        while attempt <= max_retries:
            try:
                # slice 1 item at a time for reliability
                result = molecule_client.filter().only("molecule_chembl_id")[current_index:current_index+1]
                if not result:
                    logger.info("No more results returned from API.")
                    pbar.close()
                    return all_ids

                mol = result[0]
                if "molecule_chembl_id" in mol:
                    all_ids.append(mol["molecule_chembl_id"])
                current_index += 1
                pbar.update(1)

                if resume and current_index % 100 == 0:
                    with open(checkpoint_file, "w") as f:
                        json.dump({"last_index": current_index, "collected_ids": all_ids}, f)

                time.sleep(0.05)
                break  # success, exit retry loop

            except Exception as e:
                attempt += 1
                if attempt > max_retries:
                    logger.warning(f"Max retries reached for index {current_index}. Skipping...")
                    current_index += 1
                    break
                else:
                    wait = retry_wait * (2 ** (attempt - 1)) + random.uniform(0, 0.5)
                    logger.warning(f"Error fetching molecule at index {current_index}: {e}. "
                                   f"Retry {attempt}/{max_retries} in {wait:.1f}s.")
                    time.sleep(wait)

    pbar.close()

    if resume:
        with open(checkpoint_file, "w") as f:
            json.dump({"last_index": current_index, "collected_ids": all_ids}, f)

    logger.info(f"✅ Total molecule IDs collected: {len(all_ids)}")

    if resume and cleanup and os.path.exists(checkpoint_file):
        try:
            os.remove(checkpoint_file)
            logger.info("🧹 Checkpoint file deleted after successful fetch.")
        except Exception as e:
            logger.warning(f"Failed to delete checkpoint file: {e}")

    return all_ids


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
        config={**config, "molecule_batches": [(i, batch) for i, batch in enumerate(batches)]},  # <- inject here
        input_key="molecule_batches",
        output_key="batch_index",
        filename_pattern="batch_{batch_index}.csv",
        combined_filename=config.get("output", {}).get("filename", "all_molecules_combined.csv"),
        output_dir=config.get("output", {}).get("directory", "outputs/chembl_parallel"),
        input_is_pair=True,
        use_multiprocessing=True
    )


    # Format inputs as (index, batch) tuples
    runner.inputs = [(i, batch) for i, batch in enumerate(batches)]
    return runner.run()

def _single_batch_runner(input_pair):
    batch_index, chembl_ids = input_pair
    task = get_task("retrieve_compound_metadata_batch")
    if not task:
        raise RuntimeError("Task 'retrieve_compound_metadata_batch' not registered.")
    local_config = {"molecule_chembl_ids": chembl_ids}
    return task(local_config)
