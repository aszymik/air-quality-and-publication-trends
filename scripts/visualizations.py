from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib import patches
import seaborn as sns
import pandas as pd
import geopandas as gpd

def plot_trends_for_chosen_cities(df: pd.DataFrame, chosen_years: list, chosen_cities: list):
    """Rysuje wykres trendów średnich miesięcznych PM2.5 dla wybranych dwóch miast."""

    assert len(chosen_cities) == 2, "Funkcja obsługuje dokładnie dwa miasta."
    assert len(chosen_years) == 2, "Funkcja obsługuje dokładnie dwa lata."

    plt.figure(figsize=(12, 6))

    # Wizualne rozdzielenie na trendy danych
    sns.set_style("whitegrid")

    sns.lineplot(
        data=df,
        x='Miesiąc',
        y='PM2.5',
        hue='Miejscowość',
        style='Rok',
        markers=True,
        dashes=True,
        palette={chosen_cities[0]: 'navy', chosen_cities[1]: 'darkred'},
        linewidth=2.5
    )

    plt.title(f'Trend średnich miesięcznych stężeń PM2.5: {chosen_cities[0]} vs {chosen_cities[1]} ({chosen_years[0]} vs {chosen_years[1]})', fontsize=16)
    plt.xlabel('Miesiąc', fontsize=12)
    plt.ylabel('Średnie stężenie PM2.5 [µg/m³]', fontsize=12)
    plt.xticks(range(1, 13))  
    plt.axhline(y=15, color='gray', linestyle='--', alpha=0.7, label='Norma dobowa WHO (15 µg/m³)')

    plt.legend(title='Miejscowość / Rok')
    plt.tight_layout()
    plt.show()

def plot_heatmaps_for_cities(df_means: pd.DataFrame):
    """Rysuje heatmapy miesięcznych średnich PM2.5 dla miast."""

    cities = [col for col in df_means.columns if col not in ['Data', 'Rok', 'Miesiąc']]

    n = len(cities)
    _, axes = plt.subplots(n//2, 2, figsize=(14, 2*n))
    axes = axes.flatten()

    for ax, city in zip(axes, cities):
        df_city = df_means[['Rok','Miesiąc', city]]

        # Tabela, w ktorej wiersze to lata, a kolumny miesiące
        pivot = df_city.pivot_table(
            index='Rok', 
            columns='Miesiąc', 
            values=city,
            aggfunc='mean'
        )

        sns.heatmap(pivot, annot=True, fmt=".1f", ax=ax, vmin=5, vmax=75)
        ax.set_title(f"Średnie miesięczne PM2.5 [µg/m³] – {city}")
        ax.set_ylabel("Rok")
        ax.set_xlabel("Miesiąc")
        ax.set_xticklabels(pivot.columns)

    plt.tight_layout()
    plt.show()

def plot_who_exceeding_days(selected_stations: pd.DataFrame):
    """Rysuje wykres słupkowy liczby dni przekroczeń normy WHO dla wybranych stacji."""
    plot_df = selected_stations[[2015, 2018, 2021, 2024]]

    plot_df.plot(
        kind='bar',
        figsize=(12, 7),
        width=0.4,          
        edgecolor='black',
        rot=0
    )

    plt.title('Liczba dni z przekroczeniem normy WHO (15 µg/m³) w roku\n(3 najlepsze i 3 najgorsze stacje w 2024)', fontsize=14)
    plt.xlabel('Stacja')
    plt.ylabel('Liczba dni z przekroczeniem')

    plt.axhline(y=365, color='red', linestyle=':', label='Pełny rok')
    plt.legend(title='Rok')

    plt.tight_layout()
    plt.show()

def plot_voivodeship_exceeding_days_map(
    voivodeship_counts: pd.DataFrame,
    geojson_path: str,
    years: list,
    cmap: str = 'OrRd'
):
    """Rysuje mapy województw z liczbą dni przekroczeń normy PM2.5 dla wybranych lat."""
    gdf = gpd.read_file(geojson_path)
    gdf['voiv_norm'] = gdf['nazwa'].astype(str).str.strip().str.lower()

    counts = voivodeship_counts.copy()
    counts.index = counts.index.map(lambda x: str(x).strip().lower())
    counts.index.name = 'voiv_norm'
    counts = counts.reset_index()

    vmin = counts[years].min().min()
    vmax = counts[years].max().max()

    fig, axes = plt.subplots(len(years), 1, figsize=(7, 5 * len(years)))
    if len(years) == 1:
        axes = [axes]

    for ax, year in zip(axes, years):
        data = gdf.merge(
            counts[['voiv_norm', year]],
            on='voiv_norm',
            how='left'
        )

        missing_color = 'dimgray'
        data.plot(
            column=year,
            ax=ax,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            linewidth=0.4,
            edgecolor='gray',
            missing_kwds={
                'color': missing_color,
                'edgecolor': 'gray',
                'label': 'Brak danych'
            }
        )

        data['rep_point'] = data.geometry.representative_point()
        for _, row in data.iterrows():
            value = row[year]
            if pd.notna(value):
                ax.annotate(
                    int(value),
                    xy=(row['rep_point'].x, row['rep_point'].y),
                    ha='center',
                    va='center',
                    fontsize=12
                )

        ax.set_title(f'Województwa – liczba dni > próg ({year})')
        ax.axis('off')
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=colors.Normalize(vmin=vmin, vmax=vmax))
        sm._A = []
        cb = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02, label='Liczba dni')
        cb.ax.legend(
            handles=[patches.Patch(facecolor=missing_color, edgecolor='gray', label='Brak danych')],
            loc='upper center',
            bbox_to_anchor=(0.5, -0.08),
            frameon=False
        )

    plt.tight_layout()
    plt.show()
    
