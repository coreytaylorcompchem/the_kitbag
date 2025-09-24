import requests

URL = "https://www.ebi.ac.uk/chembl/api/data/target.json?limit=1"

try:
    print("Querying ChEMBL API directly...")
    r = requests.get(URL, timeout=10)
    r.raise_for_status()
    print("✅ Response received:")
    print(r.json())
except requests.exceptions.RequestException as e:
    print(f"❌ Request failed: {e}")
