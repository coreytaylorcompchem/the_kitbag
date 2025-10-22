import requests

url = "https://bindingdb.org/rest/getLigandsByUniprot?uniprot=P00533&response=application/json"
resp = requests.get(url)
print("Status:", resp.status_code)
print("Body:", resp.text[:500], repr(resp.text))

