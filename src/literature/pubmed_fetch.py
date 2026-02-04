import argparse
import yaml
import pandas as pd
from Bio import Entrez
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description="Fetch PubMed papers based on config.")
    parser.add_argument("--year", type=int, required=True, help="Year to fetch data for")
    parser.add_argument("--config", type=str, default="config/task4.yaml", help="Path to YAML config")
    return parser.parse_args()

def fetch_pubmed_data(year, params):
    Entrez.email = params["entrez_email"]
    
    journals_q = " OR ".join([f'"{j}"[Journal]' for j in params["journals"]])
    keywords_q = " OR ".join([f'"{k}"' for k in params["keywords"]])
    query = f"({journals_q}) AND ({keywords_q}) AND {year}[pdat]"
    
    print(f"Searching for: {query}")
    
    # Pobranie ID
    with Entrez.esearch(db="pubmed", term=query, retmax=params["max_records"]) as handle:
        search_results = Entrez.read(handle)
    
    ids = search_results.get("IdList", [])
    if not ids:
        print(f"No records found for year {year}")
        return pd.DataFrame()

    # Pobranie szczegółów
    with Entrez.efetch(db="pubmed", id=",".join(ids), rettype="medline", retmode="xml") as handle:
        records = Entrez.read(handle)
    
    # Parsowanie wyników do listy słowników
    papers = []
    for article in records.get("PubmedArticle", []):
        medal = article["MedlineCitation"]["Article"]
        papers.append({
            "pmid": str(article["MedlineCitation"]["PMID"]),
            "title": medal.get("ArticleTitle", ""),
            "journal": medal.get("Journal", {}).get("Title", ""),
            "year": year,
            "doi": next((str(id) for id in article["PubmedData"]["ArticleIdList"] if id.attributes.get("IdType") == "doi"), None)
        })
    
    return pd.DataFrame(papers)

def main():
    args = parse_args()
    
    with open(args.config, 'r') as f:
        params = yaml.safe_load(f)

    output_dir = Path(f"results/literature/{args.year}")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Pobierz dane
    df = fetch_pubmed_data(args.year, params)
    
    if df.empty:
        return

    # Zapisz do CSV
    csv_path = output_dir / "pubmed_papers.csv"
    df.to_csv(csv_path, index=False)
    print(f"Saved records to {csv_path}")

    # Podsumowanie wyników
    print("\n--- Aggregations ---")

    top_journals = df['journal'].value_counts().head(5)
    print("Top Journals:\n", top_journals)
    
    summary = {
        "year": args.year,
        "total_records": len(df),
        "unique_journals": df['journal'].nunique()
    }
    print("Summary by year:", summary)


    with open(output_dir / "summary.txt", "w") as f:
        f.write(f"Summary for {args.year}:\n")
        f.write(f"Total papers: {len(df)}\n")
        f.write(f"Top journals:\n{top_journals.to_string()}\n")

if __name__ == "__main__":
    main()