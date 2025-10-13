import requests
url = "https://rest.uniprot.org/uniprotkb/search"
params = {
    "query": "interpro:IPR000276",
    "fields": "accession",
    "format": "json",
    "size": 500
}
r = requests.get(url, params=params)
print(r.status_code, r.text[:500])
