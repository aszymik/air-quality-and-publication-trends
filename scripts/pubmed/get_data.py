from Bio import Entrez
import yaml

CONFIG_PATH = "config/task4.yaml"

with open(CONFIG_PATH) as stream:
    params = yaml.safe_load(stream)

years = params["years"]
Entrez.email = params["entrez_email"]

for year in years:
    query = (
        "(" +
        " OR".join([f'"{journal}"[Journal]' for journal in params["journals"]]) +
        ") AND (" +
        " OR".join([f'"{keyword}"' for keyword in params["keywords"]]) +
        f") AND {year}[pdat]"
    )
    stream = Entrez.esearch(db="pubmed", term=query, retmax=params["max_records"])
    record = Entrez.read(stream)
    stream.close()
    print(record["Count"])
    print(len(record["IdList"]))