from pipeline.task_registry import register_task
from pipeline.logger import setup_logger
import requests

from tqdm import tqdm

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/target"
# CLASS_A_INTERPRO_IDS = {'IPR000276', 'IPR017452'}

@register_task("retrieve_protein_class_target_list", description="Retrieve UniProt IDs for protein target class.")
def retrieve_protein_class_target_list(config, data=None):
    target_class = config.get("protein_class_keyword")
    species = config.get("species")

    all_targets = []
    offset = 0
    limit = 500

    logger.info(f"Searching for targets matching keyword '{target_class}' and species: {species}")

    # Initial call to get total count
    params = {
        'pref_name__icontains': 'receptor',
        'limit': limit,
        'offset': 0
    }
    response = requests.get(BASE_URL, params=params, headers={'Accept': 'application/json'})
    response.raise_for_status()
    data = response.json()
    total_count = data['page_meta']['total_count']
    logger.info(f"Total targets found: {total_count}")

    # Process in chunks
    all_targets = data.get('targets', [])
    offset = len(all_targets)

    with tqdm(total=total_count, desc="Fetching targets", unit="targets") as pbar:
        pbar.update(len(all_targets))
        while offset < total_count:
            params['offset'] = offset
            response = requests.get(BASE_URL, params=params, headers={'Accept': 'application/json'})
            response.raise_for_status()
            data = response.json()
            targets = data.get('targets', [])
            if not targets:
                break
            all_targets.extend(targets)
            offset += len(targets)
            pbar.update(len(targets))


    logger.info(f"Retrieved {len(all_targets)} targets matching keyword '{target_class}'.")

    interpro_accession = config.get("interpro_accession") # using InterPro ID to narrow uniprot search
    logger.info(f"Using InterPro accession {interpro_accession} to assist search")

    _uniprots = []

    for target in all_targets:
        if target.get("organism") != species:
            continue

        for comp in target.get("target_components", []):
            for xref in comp.get("target_component_xrefs", []):
                if interpro_accession:
                    if xref.get("xref_src_db") == "InterPro" and xref.get("xref_id") in interpro_accession:
                        accession = comp.get("accession")
                        if accession:
                            _uniprots.append(accession)
                else:
                    accession = comp.get("accession")
                    _uniprots.append(accession)

    _uniprots = list(set(_uniprots))

    if _uniprots and interpro_accession:
        logger.info(f"Building list of {len(_uniprots)} unique UniProt IDs for targets containing '{target_class}' and '{interpro_accession}'.")
    elif _uniprots:
        logger.info(f"Building list {len(_uniprots)} unique UniProt IDs for targets containing '{target_class}'.")
    else:
        logger.warning("No UniProt IDs found for the filtered targets.")

    config["uniprot_ids"] = _uniprots

    return _uniprots
