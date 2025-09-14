import requests

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/target"

CLASS_A_INTERPRO_IDS = {'IPR000276', 'IPR017452'}

def fetch_targets_with_receptor():
    all_targets = []
    offset = 0
    limit = 500  # Reduced chunk size to 500
    while True:
        params = {
            'pref_name__icontains': 'receptor',
            'limit': limit,
            'offset': offset
        }
        print(f"Fetching targets {offset} to {offset + limit}...")
        response = requests.get(BASE_URL, params=params, headers={'Accept': 'application/json'})
        response.raise_for_status()
        data = response.json()
        targets = data.get('targets', [])
        if not targets:
            break
        all_targets.extend(targets)
        offset += len(targets)
        if offset >= data['page_meta']['total_count']:
            break
    return all_targets

def is_class_a_gpcr(target):
    for comp in target.get('target_components', []):
        for xref in comp.get('target_component_xrefs', []):
            if xref.get('xref_src_db') == 'InterPro' and xref.get('xref_id') in CLASS_A_INTERPRO_IDS:
                return True
    return False

def extract_uniprot_ids(target):
    uniprot_ids = []
    for comp in target.get('target_components', []):
        accession = comp.get('accession')
        if accession:
            uniprot_ids.append(accession)
    return uniprot_ids

def main():
    targets = fetch_targets_with_receptor()
    print(f"Total targets with 'receptor' in name: {len(targets)}")
    
    class_a_gpcrs = []
    for target in targets:
        if is_class_a_gpcr(target):
            class_a_gpcrs.append({
                'target_chembl_id': target.get('target_chembl_id'),
                'pref_name': target.get('pref_name'),
                'uniprot_accessions': extract_uniprot_ids(target)
            })
    
    print(f"Filtered Class A GPCR targets count: {len(class_a_gpcrs)}")
    for entry in class_a_gpcrs:
        print(entry['target_chembl_id'], entry['pref_name'], entry['uniprot_accessions'])

if __name__ == "__main__":
    main()
