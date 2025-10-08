from pipeline.task_registry import register_task
from pipeline.logger import setup_logger
import requests

from tqdm import tqdm

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/target"
# CLASS_A_INTERPRO_IDS = {'IPR000276', 'IPR017452'}

@register_task("retrieve_protein_class_target_list", 
               category='Bioactivity',
               description="Retrieve UniProt IDs for protein target class.")
def retrieve_protein_class_target_list(config, data=None):
    strict_interpro = config.get("interpro_protein_class")
    flexible_keyword = config.get("protein_class_keyword")
    interpro_accession = config.get("interpro_accession_refinement")
    species = config.get("species")

        # --- Sanity checks ---
    if strict_interpro:
        if flexible_keyword or interpro_accession:
            raise ValueError(
                "You have specified 'protein_class' for strict search, "
                "but also provided 'protein_class_keyword' or 'interpro_accession'. "
                "These are only valid for flexible search. Please choose one mode."
            )
    elif not flexible_keyword:
        raise ValueError(
            "You must specify either 'protein_class' (strict) "
            "or 'protein_class_keyword' (flexible) in your YAML."
        )
    _uniprots = []

        # STRICT SEARCH MODE
    if strict_interpro:
        logger.info(f"Running strict search for InterPro accession: {strict_interpro}")

        offset = 0
        limit = 500
        all_targets = []

        params = {
            'limit': limit,
            'offset': offset
        }

        response = requests.get(BASE_URL, params=params, headers={'Accept': 'application/json'})
        response.raise_for_status()
        data = response.json()
        total_count = data['page_meta']['total_count']
        logger.info(f"Total targets found in ChEMBL: {total_count}")

        all_targets.extend(data.get('targets', []))
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

        logger.info(f"Filtering targets for species '{species}' and InterPro accession '{strict_interpro}'")
        for target in all_targets:
            if target.get("organism") != species:
                continue

            for comp in target.get("target_components", []):
                for xref in comp.get("target_component_xrefs", []):
                    if xref.get("xref_src_db") == "InterPro" and xref.get("xref_id") == strict_interpro:
                        accession = comp.get("accession")
                        if accession:
                            _uniprots.append(accession)

        logger.info(f"Strict search retrieved {len(_uniprots)} UniProt IDs for InterPro: {strict_interpro}")

    # FLEXIBLE SEARCH MODE 
    else:
        if interpro_accession:
            logger.info(f"Running flexible search for keyword '{flexible_keyword}', InterPro accession '{interpro_accession}', and species: {species}")
        else:
            logger.info(f"Running flexible search for keyword '{flexible_keyword}' and species: {species}")


        offset = 0
        limit = 500
        params = {
            'pref_name__icontains': flexible_keyword,
            'limit': limit,
            'offset': 0
        }

        response = requests.get(BASE_URL, params=params, headers={'Accept': 'application/json'})
        response.raise_for_status()
        data = response.json()
        total_count = data['page_meta']['total_count']
        logger.info(f"Total targets found: {total_count}")

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

        logger.info(f"Filtering {len(all_targets)} targets with optional InterPro filter: {interpro_accession}")

        for target in all_targets:
            if target.get("organism") != species:
                continue

            for comp in target.get("target_components", []):
                accession = comp.get("accession")
                if not accession:
                    continue

                if interpro_accession:
                    for xref in comp.get("target_component_xrefs", []):
                        if xref.get("xref_src_db") == "InterPro" and xref.get("xref_id") in interpro_accession:
                            _uniprots.append(accession)
                            break
                else:
                    _uniprots.append(accession)

        logger.info(f"Flexible search retrieved {len(_uniprots)} UniProt IDs using keyword '{flexible_keyword}'.")

    # Remove duplicates
    _uniprots = list(set(_uniprots))

    if not _uniprots:
        logger.warning("No UniProt IDs found.")
    else:
        logger.info(f"Final UniProt ID count: {len(_uniprots)}")

    config["uniprot_ids"] = _uniprots
    return _uniprots

