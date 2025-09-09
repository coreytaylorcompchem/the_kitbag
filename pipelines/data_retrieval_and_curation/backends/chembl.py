import requests
from .base import DataBackend

class ChemblBackend(DataBackend):
    def fetch_bioactivity(self, compound_id):
        url = f"https://www.ebi.ac.uk/chembl/api/data/activity.json?molecule_chembl_id={compound_id}"
        response = requests.get(url)
        return response.json().get('activities', [])

    def fetch_structure(self, compound_id):
        url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{compound_id}.json"
        response = requests.get(url)
        return response.json()
