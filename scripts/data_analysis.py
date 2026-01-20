import pandas as pd

def get_monthly_means_for_stations(df: pd.DataFrame) -> pd.DataFrame:
    """Oblicza miesięczne średnie wartości PM2.5 dla każdej stacji.
    Arguments:
        df: DataFrame z danymi PM2.5, gdzie kolumny to kody stacji, a indeks to daty.
    Returns:
        DataFrame z miesięcznymi średnimi wartościami PM2.5 dla każdej stacji.
    """
    # Konwersja tylko kolumn stacji
    df_num = df.drop(columns=[('Data', '')]).apply(
        pd.to_numeric, errors="coerce"
    )
    monthly_means = (
        df_num
        .groupby([
            df[('Data', '')].dt.year.rename("Rok"),
            df[('Data', '')].dt.month.rename("Miesiąc"),
        ])
        .mean()
    )
    return monthly_means

def get_chosen_monthly_means(df: pd.DataFrame, chosen_years: list, chosen_cities: list) -> pd.DataFrame:
    """Oblicza miesięczne średnie wartości PM2.5 dla wybranych miast i lat.
    Arguments:
        df: DataFrame z danymi PM2.5, gdzie kolumny to (miejscowość, kod stacji), a indeks to daty.
        chosen_years: lista lat do uwzględnienia.
        chosen_cities: lista miast do uwzględnienia.
    Returns:
        DataFrame z miesięcznymi średnimi wartościami PM2.5 dla wybranych miast i lat.
    """
    # Filtrowanie po miastach i latach
    city_cols = [col for col in df.columns if col[0] in chosen_cities]
    year_rows = df[('Data', '')].dt.year.isin(chosen_years)
    
    # Tworzymy format długi do grupowania
    df_long = df.loc[year_rows, [('Data', '')] + city_cols].melt(
        id_vars=[('Data', '')],
        var_name='Miejscowość',
        value_name='PM2.5'
    )

    # Średnie miesięczne po miastach
    df_monthly = (
        df_long.groupby([
            df_long[('Data', '')].dt.year.rename('Rok'),
            df_long[('Data', '')].dt.month.rename('Miesiąc'),
            'Miejscowość'
        ])['PM2.5']
        .mean()
        .reset_index()
    )
    return df_monthly

def get_monthly_means_for_cities(df: pd.DataFrame) -> pd.DataFrame:
    """Oblicza miesięczne średnie PM2.5 uśrednione dla wszystkich stacji miasta.
    Arguments:
        df: DataFrame z danymi PM2.5, gdzie kolumny to (miejscowość, kod stacji), a indeks to daty.
    Returns:
        DataFrame z miesięcznymi średnimi wartościami PM2.5 uśrednionymi dla każdego miasta.
    """
    # Wybieramy kolumny z danymi
    data_cols = [col for col in df.columns if col[0] != 'Data']
    df_num = df[data_cols].apply(pd.to_numeric, errors='coerce')

    # Grupowanie po miejscowości (poziom 0 indeksu) i uśrednienie
    df_city = df_num.T.groupby(level=0).mean().T

    # Uśrednienie po miesiącu w każdym roku
    df_means = (
        df_city.groupby([
            df[('Data', '')].dt.year.rename('Rok'),
            df[('Data', '')].dt.month.rename('Miesiąc')
        ])
        .mean()
        .reset_index()
    )
    return df_means

def get_who_norm_exceeding_days(df: pd.DataFrame) -> pd.DataFrame:
    """Zwraca liczbę dni w miesiącu, w których średnie dzienne PM2.5 przekroczyły normę WHO (15 µg/m³).
    Arguments:
        df: DataFrame z danymi PM2.5, gdzie kolumny to (miejscowość, kod stacji), a indeks to daty.
    Returns:
        DataFrame z liczbą dni w miesiącu przekraczających normę WHO dla każdej stacji.
    """
    # Ustawiamy indeks czasowy i resamplujemy dane do częstotliwości dziennej ('D')
    df[('Data', '')] = pd.to_datetime(df[('Data', '')], format="mixed")
    # Ustawiamy datę jako indeks
    df_daily = df.set_index(('Data', ''))
    daily_means = df_daily.resample('D').mean()

    exceeded_mask = daily_means > 15
    yearly_counts = exceeded_mask.groupby(daily_means.index.year).sum()

    return yearly_counts.T  # stacje w wierszach, lata w kolumnach

def get_max_and_min_k_stations(yearly_counts: pd.DataFrame, chosen_year: int, k: int=3) -> pd.DataFrame:
    """Zwraca DataFrame z k stacjami o najwyższych i najniższych liczbach dni 
    z przekroczeniem normy WHO w wybranym roku.
    Arguments:
        yearly_counts: DataFrame z liczbą dni przekroczeń normy WHO, gdzie wiersze to stacje, a kolumny to lata.
        chosen_year: rok do analizy.
        k: liczba stacji do zwrócenia z najwyższymi i najniższymi wartościami.
    Returns:
        DataFrame z 2k stacjami (k najwyższych i k najniższych).
    """
    if chosen_year not in yearly_counts.columns:
        return pd.DataFrame()
        
    sorted_results = yearly_counts.sort_values(by=chosen_year)
    return pd.concat([sorted_results.head(k), sorted_results.tail(k)])

def get_voivodeship_exceeding_days(df: pd.DataFrame, code_to_voivodeship: dict, threshold: float=15) -> pd.DataFrame:
    """Zwraca liczbę dni w roku, w których średnie dzienne PM2.5 przekroczyły próg
    dla dowolnej stacji w danym województwie.
    Arguments:
        df: DataFrame z danymi PM2.5, gdzie kolumny to (miejscowość, kod stacji).
        code_to_voivodeship: słownik mapujący kody stacji na województwa.
        threshold: próg stężenia PM2.5 (µg/m³).
    Returns:
        DataFrame z liczbą dni przekroczeń dla województw (wiersze) i lat (kolumny).
    """
    df = df.copy()
    df[('Data', '')] = pd.to_datetime(df[('Data', '')], format="mixed")

    station_cols = [col for col in df.columns if col != ('Data', '')]
    df_daily = df.set_index(('Data', ''))
    df_num = df_daily[station_cols].apply(pd.to_numeric, errors='coerce')
    daily_means = df_num.resample('D').mean()

    exceeded = daily_means > threshold

    if isinstance(exceeded.columns, pd.MultiIndex):
        station_codes = [col[1] for col in exceeded.columns]
    else:
        station_codes = list(exceeded.columns)

    voivodeships = [
        code_to_voivodeship.get(code, 'Nieznane')
        for code in station_codes
    ]
    exceeded.columns = pd.MultiIndex.from_arrays(
        [voivodeships, station_codes],
        names=['Województwo', 'Kod stacji']
    )

    exceeded_voiv = exceeded.T.groupby(level=0).any().T
    yearly_counts = exceeded_voiv.groupby(exceeded_voiv.index.year).sum()
    return yearly_counts.T
