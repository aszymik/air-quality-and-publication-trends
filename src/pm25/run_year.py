import argparse
import json
import os
import yaml

from load_data import (
    get_metadata,
    get_code_mappings,
    download_and_preprocess_data
)
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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--year", type=int, required=True)
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    year = args.year
    if year < 2006 or year > 2024:
        print(f'No data for the chosen year ({year}). Please select a year between 2006 and 2024.')
        return
    outdir = f"results/pm25/{year}"
    os.makedirs(outdir, exist_ok=True)

    with open(args.config) as f:
        config = yaml.safe_load(f)
    
    cities = config["cities"]
    print(cities)
    gios_id_data = json.loads(open("data/gios_ids.json").read())
    gios_ids = {int(k): v for k, v in gios_id_data.items()}

    # Pobierz metadane
    metadata = get_metadata()
    old_to_new_code, code_to_city, code_to_voivodeship = get_code_mappings(metadata)

    # Pobierz dane
    df = download_and_preprocess_data(
        year=args.year,
        gios_id=gios_ids[year],
        gios_filename=f'{year}_PM25_1g.xlsx',
        code_to_city=code_to_city,
        old_to_new_code=old_to_new_code
    )
    os.makedirs(f'data/tables/{year}', exist_ok=True)
    df.to_csv(f'data/tables/{year}/{year}_data.csv')

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

    geojson_path = 'data/wojewodztwa-min.geojson'
    plot_voivodeship_exceeding_days_map(voiv_counts, geojson_path, [year], f"{outdir}/figures")

if __name__ == "__main__":
    main()