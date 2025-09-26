import os
import requests
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio.PDB import PDBList

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query?json="

def get_pdb_ids_for_uniprot(uniprot_id):
    url = "https://search.rcsb.org/rcsbsearch/v2/query"

    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": uniprot_id
            }
        },
        "return_type": "entry"
    }

    try:
        response = requests.post(url, json=query, timeout=10)

        # Check for HTTP errors
        if response.status_code != 200:
            logger.error(f"RCSB API query failed with status {response.status_code}")
            logger.debug(f"Response text: {response.text}")
            return []

        try:
            response_json = response.json()
        except ValueError as e:
            logger.error("❌ Failed to parse JSON from RCSB API response.")
            logger.debug(f"Response content: {response.text}")
            return []

        pdb_ids = [result["identifier"] for result in response_json.get("result_set", [])]
        return pdb_ids

    except requests.RequestException as e:
        logger.error(f"❌ Exception while querying RCSB API: {e}")
        return []


def download_pdb(pdb_id, output_dir, pdbl):
    try:
        filepath = pdbl.retrieve_pdb_file(pdb_id, pdir=output_dir, file_format="pdb", overwrite=True)
        logger.debug(f"Downloaded: {pdb_id}")
        return filepath
    except Exception as e:
        logger.warning(f"Failed to download {pdb_id}: {e}")
        return None

@register_task("retrieve_pdbs", category="PDB", description="Retrieve all PDBs for a UniProt ID.")
def retrieve_pdbs(config, data=None):
    uniprot_id = config["uniprot_id"]
    output_dir = config.get("output", {}).get("directory", "outputs/pdb_downloads")
    os.makedirs(output_dir, exist_ok=True)

    pdb_ids = get_pdb_ids_for_uniprot(uniprot_id)
    if not pdb_ids:
        logger.warning(f"No PDBs found for UniProt ID {uniprot_id}")
        return []

    logger.info(f"Found {len(pdb_ids)} PDBs for UniProt ID {uniprot_id}")

    pdbl = PDBList()
    downloaded = []

    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = {
            executor.submit(download_pdb, pdb_id, output_dir, pdbl): pdb_id
            for pdb_id in pdb_ids
        }
        for future in as_completed(futures):
            result = future.result()
            if result:
                downloaded.append(result)

    logger.info(f"{len(downloaded)} PDBs downloaded successfully.")
    return downloaded