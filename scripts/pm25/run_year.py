import argparse
import yaml
import pandas as pd
from pathlib import Path

from load_data import (
    get_metadata,
    get_code_mappings,
    download_and_preprocess_data
)
from data_analysis import (
    get_monthly_means_for_cities,
    get_who_norm_exceeding_days
)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--year", type=int, required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    # ---- metadata ----
    metadata = get_metadata()
    old_to_new, code_to_city, _ = get_code_mappings(metadata)

    # ---- data ----
    df = download_and_preprocess_data(
        year=args.year,
        gios_id="1234",              # to do
        gios_filename="pm25.xlsx",   # to do
        code_to_city=code_to_city,
        old_to_new_code=old_to_new
    )

    # ---- daily / monthly ----
    daily_means = get_monthly_means_for_cities(df)
    daily_means.to_csv(outdir / "daily_means.csv", index=False)

    # ---- exceedance ----
    exceed = get_who_norm_exceeding_days(df)
    exceed.to_csv(outdir / "exceedance_days.csv")

if __name__ == "__main__":
    main()
