import argparse
import json
import os
import pandas as pd

from load_data import (
    get_code_mappings,
    download_and_preprocess_data
)

GIOS_ID_PATH = "data/gios_ids.json"

def main():
    parser = argparse.ArgumentParser(description="Pobiera dane PM2.5 z danego roku.")
    parser.add_argument(
        "--year", 
        type=int,
        required=True, 
        help="Wybrany rok (z zakresu 2006-2024)"
    )
    parser.add_argument(
        "--metadata", 
        required=True, 
        help="Ścieżka do pliku CSV z metadanymi stacji (np. data/tables/metadata.csv)"
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Katalog, w którym zostanie zapisany plik z danymi."
    )
    args = parser.parse_args()
    year = args.year
    if year < 2006 or year > 2024:
        print(f'Brak danych dla wybranego roku ({year}). Proszę wybrać rok z zakresu 2006–2024.')
        return

    metadata = pd.read_csv(args.metadata, index_col=0) 
    old_to_new_code, code_to_city, _ = get_code_mappings(metadata)
    gios_id_data = json.loads(open(GIOS_ID_PATH).read())
    gios_ids = {int(k): v for k, v in gios_id_data.items()}

    print("Pobieranie danych...")
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    data_csv_path = f'{output_dir}/{year}.csv'
    
    try:
        df = download_and_preprocess_data(
                year=args.year,
                gios_id=gios_ids[year],
                gios_filename=f'{year}_PM25_1g.xlsx',
                code_to_city=code_to_city,
                old_to_new_code=old_to_new_code
            )
    except KeyError:
        df = download_and_preprocess_data(
                year=args.year,
                gios_id=gios_ids[year],
                gios_filename=f'{year}_PM2.5_1g.xlsx',
                code_to_city=code_to_city,
                old_to_new_code=old_to_new_code
            )
    df.to_csv(data_csv_path, index=None)

if __name__ == "__main__":
    main()