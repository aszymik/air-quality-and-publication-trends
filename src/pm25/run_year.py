import argparse
import json
import os
import yaml
import pandas as pd

from load_data import get_code_mappings

from data_analysis import (
    get_monthly_means_for_stations,
    get_chosen_monthly_means,
    get_monthly_means_for_cities,
    get_who_norm_exceeding_days,
    get_max_and_min_k_stations,
    get_voivodeship_exceeding_days
)

from visualizations import (
    plot_trends_for_chosen_cities,
    plot_heatmaps_for_cities,
    plot_who_exceeding_days,
    plot_voivodeship_exceeding_days_map
)

GEOJSON_PATH = "data/wojewodztwa-min.geojson"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--year", type=int, required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--data_dir", required=True)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    year = args.year
    if year < 2006 or year > 2024:
        print(f'No data for the chosen year ({year}). Please select a year between 2006 and 2024.')
        return
    outdir = f"{args.output_dir}/{year}"
    os.makedirs(outdir, exist_ok=True)

    with open(args.config) as f:
        config = yaml.safe_load(f)
    cities = config["cities"]

    # Wczytanie danych i metadanych
    metadata = pd.read_csv(args.metadata, index_col=0)
    _, _, code_to_voivodeship = get_code_mappings(metadata)
    df = pd.read_csv(f"{args.data_dir}/{year}.csv", index_col=0)

    # Obliczanie średnich miesięcznych
    monthly_means = get_monthly_means_for_stations(df)
    monthly_means.to_csv(f"{outdir}/monthly_means_stations.csv", index=False)

    # Trendy dla wybranych miast
    year = int(year)
    df_plot = get_chosen_monthly_means(df, [year], cities)
    plot_trends_for_chosen_cities(df_plot, year, cities, f"{outdir}/figures")

    # Średnie miesięczne
    df_means = get_monthly_means_for_cities(df)
    df_means.to_csv(f"{outdir}/monthly_means_cities.csv", index=False)
    plot_heatmaps_for_cities(df_means, f"{outdir}/figures")

    # Dni z przekroczeniem normy WHO
    exceed = get_who_norm_exceeding_days(df)
    exceed.to_csv(f"{outdir}/exceedance_days.csv")
    top_stations = get_max_and_min_k_stations(exceed, chosen_year=year, k=3)
    plot_who_exceeding_days(top_stations, f"{outdir}/figures")

    # Mapa dni z przekroczeniem normy dla województw
    voiv_counts = get_voivodeship_exceeding_days(df, code_to_voivodeship, threshold=15)
    voiv_counts.to_csv(f"{outdir}/voivodeship_exceedance_days.csv")

    plot_voivodeship_exceeding_days_map(voiv_counts, GEOJSON_PATH, [year], f"{outdir}/figures")

if __name__ == "__main__":
    main()