import time
import os
import requests

import pandas as pd

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

    downloaded = []  # list of paths to downloaded .ent files

    # After downloading:
    renamed_files = rename_ent_files(output_dir)

    logger.info(f"{len(renamed_files)} PDBs downloaded and renamed successfully.")
    return renamed_files
    return downloaded

# --- New task to analyze InterPro family and get UniProt IDs + stats ---

def rename_ent_files(directory):
    renamed_files = []
    for filename in os.listdir(directory):
        if filename.endswith(".ent") and filename.startswith("pdb"):
            old_path = os.path.join(directory, filename)
            new_name = filename[3:-4] + ".pdb"  # Strip "pdb" prefix and change extension
            new_path = os.path.join(directory, new_name)
            os.rename(old_path, new_path)
            logger.info(f"Renamed {filename} -> {new_name}")
            renamed_files.append(new_path)
    return renamed_files

def get_uniprots_for_interpro(interpro_id, max_pages=200):
    """
    Retrieve UniProt accessions annotated with a given InterPro ID (e.g. IPR000276)
    using the UniProt REST API with the correct cross-reference query syntax.
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"xref:InterPro-{interpro_id}",
        "fields": "accession",
        "format": "json",
        "size": 500
    }

    uniprots = set()
    page = 1
    url = base_url

    logger.debug(f"Fetching UniProt accessions for InterPro {interpro_id} from UniProt API...")

    while url and page <= max_pages:
        try:
            if page == 1:
                r = requests.get(url, params=params, timeout=60)
            else:
                r = requests.get(url, timeout=60)  # no params here, next_url already contains them

            r.raise_for_status()
            data = r.json()
        except requests.RequestException as e:
            logger.error(f"❌ Failed fetching page {page} for {interpro_id}: {e}")
            break

        # Extract UniProt accessions
        for result in data.get("results", []):
            acc = result.get("primaryAccession")
            if acc:
                uniprots.add(acc)

        # Handle pagination
        next_url = None
        for link in data.get("links", []):
            if link.get("rel") == "next":
                next_url = link.get("href")
                break

        url = next_url
        params = None
        page += 1

        if page % 10 == 0:
            logger.info(f"Processed {page} pages ({len(uniprots)} UniProt IDs so far)")
        time.sleep(0.2)  # be polite

    logger.debug(f"Retrieved {len(uniprots)} UniProt accessions for {interpro_id}")
    return sorted(uniprots)



def get_pdb_resolution_stats(pdb_ids):
    """
    Query resolution info for a list of PDB IDs and compute mean resolution.
    Only includes structures with available resolution.
    """
    if not pdb_ids:
        return None

    resolutions = []
    for pdb_id in pdb_ids:
        try:
            url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
            resp = requests.get(url, timeout=10)
            if resp.status_code != 200:
                continue
            data = resp.json()
            reso = data.get("rcsb_entry_info", {}).get("resolution_combined")
            if reso:
                # resolution_combined is a list; take the first if available
                if isinstance(reso, list):
                    resolutions.append(reso[0])
                else:
                    resolutions.append(reso)
        except Exception as e:
            logger.debug(f"Failed to fetch resolution for {pdb_id}: {e}")

    if resolutions:
        mean_res = sum(resolutions) / len(resolutions)
        return mean_res
    else:
        return None


@register_task("analyse_pdbs_by_interpro_accession", category="PDB", description='Check how many PDBs are available for a given InterPro ID.')
def analyse_pdbs_by_interpro_accession(config, data=None):
    interpro_ids = config.get("interpro_id", [])
    
    # Normalize to list
    if isinstance(interpro_ids, str):
        interpro_ids = [interpro_ids]
    elif not isinstance(interpro_ids, (list, tuple)):
        raise ValueError(f"Unexpected type for interpro_id: {type(interpro_ids)}")
    
    all_uniprots = []
    for ipr in interpro_ids:
        logger.debug(f"Fetching UniProt accessions for InterPro {ipr} from InterPro API...")
        ids = get_uniprots_for_interpro(ipr)
        all_uniprots.extend(ids)

    if not all_uniprots:
        logger.warning("No UniProt entries found for InterPro analysis.")
        return {"uniprots": []}

    logger.debug(f"Total {len(all_uniprots)} UniProt IDs collected across {len(interpro_ids)} InterPro entries.")
    return {"uniprots": all_uniprots}
