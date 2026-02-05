import argparse
import os
import pandas as pd
from load_data import get_metadata

def main():
    parser = argparse.ArgumentParser(description="Pobiera metadane stacji.")
    parser.add_argument(
        "--output", 
        required=True, 
        help="Ścieżka do pliku wyjściowego CSV (np. data/tables/metadata.csv)"
    )
    args = parser.parse_args()

    print("Pobieranie metadanych...")
    metadata = get_metadata()

    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    metadata.to_csv(args.output)
    print(f"Metadane zapisane do {args.output}")

if __name__ == '__main__':
    main()