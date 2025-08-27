import requests

url = "https://search.rcsb.org/rcsbsearch/v2/query?json="
# Corpo da query em JSON, filtrando por organismo Homo sapiens
query = {
    "query": {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name",
            "operator": "exact_match",
            "value": "Homo sapiens"
        }
    },
    "return_type": "entry",
    "request_options": {
        "paginate": {
            "start": 0,
            "rows": 100
        }
    },
}

r = requests.post(url, json=query)

print(r.status_code)

data = r.json()
codes = [entry['identifier'] for entry in data.get('result_set', [])]

with open("first_100_homo_sapiens.txt", "w") as f:
    for c in codes:
        f.write(c)
        f.write("\n")
