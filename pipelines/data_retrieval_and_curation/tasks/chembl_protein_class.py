from pipeline.task_registry import register_task
from pipeline.logger import setup_logger
import requests

from tqdm import tqdm

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/target"
CLASS_A_INTERPRO_IDS = {'IPR000276', 'IPR017452'}

@register_task("retrieve_protein_class_target_list", description="Retrieve UniProt IDs for Class A GPCRs.")
def retrieve_protein_class_target_list(config, data=None):
    species = config.get("species", "Homo sapiens")

    all_targets = []
    offset = 0
    limit = 500

    logger.info(f"🔍 Searching for targets with 'receptor' in name for species: {species}")

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
    logger.info(f"📊 Total receptor-like targets found: {total_count}")

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


    logger.info(f"📦 Retrieved {len(all_targets)} targets containing 'receptor' in pref_name.")

    class_a_gpcr_uniprots = []

    for target in all_targets:
        if target.get("organism") != species:
            continue

        for comp in target.get("target_components", []):
            for xref in comp.get("target_component_xrefs", []):
                if xref.get("xref_src_db") == "InterPro" and xref.get("xref_id") in CLASS_A_INTERPRO_IDS:
                    accession = comp.get("accession")
                    if accession:
                        class_a_gpcr_uniprots.append(accession)

    class_a_gpcr_uniprots = list(set(class_a_gpcr_uniprots))  # remove duplicates

    if class_a_gpcr_uniprots:
        logger.info(f"✅ Found {len(class_a_gpcr_uniprots)} UniProt IDs for Class A GPCRs.")
    else:
        logger.warning("⚠ No UniProt IDs found for the filtered Class A GPCR targets.")

    # Update config so ParallelRunner can access the IDs
    config["uniprot_ids"] = class_a_gpcr_uniprots

    return class_a_gpcr_uniprots
