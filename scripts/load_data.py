import pandas as pd
import requests
import zipfile
import io
from typing import Tuple


def download_gios_archive(year, gios_id, filename):
    """Pobiera podane archiwum.
    Arguments:
        year: rok danych (używane tylko do komunikatów o błędach).
        gios_id: ID archiwum w serwisie GIOŚ.
        filename: nazwa pliku wewnątrz archiwum ZIP.
    Returns:
        DataFrame z danymi PM2.5."""
    # Pobranie archiwum ZIP do pamięci
    gios_archive_url = "https://powietrze.gios.gov.pl/pjp/archives/downloadFile/"
    url = f"{gios_archive_url}{gios_id}"
    response = requests.get(url)
    response.raise_for_status()  # jeśli błąd HTTP, zatrzymaj
    
    # Otwórz zip w pamięci
    with zipfile.ZipFile(io.BytesIO(response.content)) as z:
        # znajdź właściwy plik z PM2.5
        if not filename:
            print(f"Błąd: nie znaleziono {filename}.")
        else:
            # wczytaj plik do pandas
            with z.open(filename) as f:
                try:
                    df = pd.read_excel(f, header=None)
                except Exception as e:
                    print(f"Błąd przy wczytywaniu {year}: {e}")
    return df

def download_metadata():
    """Pobiera metadane o stacjach.
    Returns:
        DataFrame z metadanymi o stacjach."""
    url = "https://powietrze.gios.gov.pl/pjp/archives/downloadFile/622"
    response = requests.get(url)
    response.raise_for_status()
    df = pd.read_excel(io.BytesIO(response.content), header=0)
    return df

def get_metadata():
    """Pobiera i preprocesuje metadane o stacjach.
    Returns:
        DataFrame z metadanymi o stacjach."""
    metadata = download_metadata()
    metadata['Stary kod stacji'] = metadata['Stary Kod stacji \n(o ile inny od aktualnego)']
    metadata = metadata[['Kod stacji', 'Stary kod stacji', 'Miejscowość', 'Województwo']]
    metadata['Stary kod stacji'] = metadata['Stary kod stacji'].replace({' ': pd.NA})  # zamieniamy spacje na nan
    return metadata

def get_code_mappings(metadata: pd.DataFrame) -> Tuple[dict, dict, dict]:
    """Tworzy słowniki mapujące stare kody na nowe i kody na miasta.
    Arguments:
        metadata: DataFrame z metadanymi o stacjach.
    Returns:
        Tuple zawierający dwa słowniki:
        - old_to_new_code: mapowanie starych kodów na nowe kody stacji.
        - code_to_city: mapowanie kodów stacji na nazwy miast.
        - code_to_voivodeship: mapowanie kodów stacji na nazwy województw."""
    metadata = metadata.copy()

    old_to_new_code = {}  # słownik mapujący stare kody na nowe

    for _, row in metadata.dropna(subset=['Stary kod stacji']).iterrows():
        old_codes = str(row['Stary kod stacji']).split(',')  # rozdzielamy stare kody po przecinku
        for code in old_codes:
            code = code.strip()
            if code:
                old_to_new_code[code] = row['Kod stacji']
    
    code_to_city = ( # słownik mapujący kody na miasta
        metadata.set_index('Kod stacji')['Miejscowość']
        .to_dict()
    )
    code_to_voivodeship = (
        metadata.set_index('Kod stacji')['Województwo']
        .to_dict()
    )
    return old_to_new_code, code_to_city, code_to_voivodeship

def rename_columns(df: pd.DataFrame, old_to_new_code: dict) -> pd.DataFrame:
    """Zmienia nazwy kolumn na nowe kody stacji.
    Arguments:
        df: DataFrame z danymi PM2.5.
        old_to_new_code: słownik mapujący stare kody na nowe kody stacji.
    Returns:
        DataFrame z zaktualizowanymi nazwami kolumn."""
    new_col_names = []
    for col in df.columns:
        if col in old_to_new_code:
            new_col_names.append(old_to_new_code[col])
        else:
            new_col_names.append(col)
    df.columns = new_col_names
    return df

def add_multiindex(df, code_dict: dict) -> pd.DataFrame:
    """Dodaje miasto do multiindeksu.
    Arguments:
        df: DataFrame z danymi PM2.5, gdzie kolumny to kody stacji.
        code_dict: słownik mapujący kody stacji na nazwy miast.
    Returns:
        DataFrame z MultiIndexem (miejscowość, kod stacji)."""
    # Pomijamy kolumnę 'Data'
    data_col = df['Data']
    data_values = df.drop(columns=['Data'])
    
    # Tworzymy MultiIndex na podstawie słownika metadanych
    for col in data_values.columns:
        if col not in code_dict:  # check, czy wszystkie kody są w słowniku
            print(col)

    tuples = [(code_dict.get(col, 'Nieznane'), col) for col in data_values.columns]
    multi_index = pd.MultiIndex.from_tuples(tuples, names=['Miejscowość', 'Kod stacji'])
    data_values.columns = multi_index
    
    # Dodajemy z powrotem kolumnę Data
    data_values.insert(0, 'Data', data_col)
    return data_values

def change_midnight_measurements(df: pd.DataFrame) -> pd.DataFrame:
    """Przesuwa pomiary o północy o jeden dzień wstecz.
    Arguments:
        df: DataFrame z danymi PM2.5, gdzie jedna z kolumn to 'Data'.
    Returns:
        DataFrame z przesuniętymi pomiarami o północy."""
    # df['Data'] = pd.to_datetime(df['Data'], format='%m/%d/%y %H:%M', errors='coerce')
    df['Data'] = pd.to_datetime(df['Data'], errors='coerce')
    df['Data'] = df['Data'].dt.floor('min')  # usuwamy mikrosekundy

    # Przesuwamy pomiary o północy o -1
    midnight_dates = df['Data'].dt.time == pd.Timestamp('00:00:00').time()
    df.loc[midnight_dates, 'Data'] = df.loc[midnight_dates, 'Data'] - pd.Timedelta(seconds=1)

    return df

def download_and_preprocess_data(year: int, gios_id: str, gios_filename: str, code_to_city: dict, old_to_new_code: dict, header_index: int=0) -> pd.DataFrame:
    """Pobiera i przygotowuje dane z archiwum GIOŚ dla podanego roku.
    Arguments:
        year: rok danych.
        gios_id: ID archiwum w serwisie GIOŚ.
        gios_filename: nazwa pliku wewnątrz archiwum ZIP.
        code_to_city: słownik mapujący kody stacji na nazwy miast.
        old_to_new_code: słownik mapujący stare kody na nowe kody stacji.
        header_index: indeks wiersza nagłówka w pliku Excel.
    Returns:
        DataFrame z przetworzonymi danymi PM2.5."""
    df = download_gios_archive(year, gios_id, gios_filename)

    # Zmieniamy nazwy kolumn na kody stacji
    col_names = df.iloc[header_index]
    col_names[0] = 'Data'
    df.columns = col_names
    df = df.drop(df.index[header_index])  # usuwamy wiersz nagłówka

    # Ujednolicamy nazwy kolumn
    df = rename_columns(df, old_to_new_code)

    # Usuwamy niepotrzebne wiersze nagłówkowe
    df = df[~df['Data'].isin(('Nr', 'Wskaźnik', 'Czas uśredniania', 'Jednostka', 'Kod stanowiska', 'Czas pomiaru'))]

    # Dla 2018 roku (inny format) konwertujemy dane na numeryczne
    if year == 2018:
        value_cols = df.columns.difference(['Data'])
        df[value_cols] = (
            df[value_cols]
            .astype(str)
            .apply(lambda col: col.str.replace(",", ".", regex=False))
            .astype(float)
        )

    # Przesuwamy pomiary z północy
    df = change_midnight_measurements(df)

    # Dodajemy multiindex: (miejscowość, kod stacji)
    df = add_multiindex(df, code_to_city)
    return df

def join_data_on_common_stations(dfs: list[pd.DataFrame]) -> tuple[pd.DataFrame, list]:
    """Łączy DataFrame'y, zachowując tylko wspólne stacje.
    Arguments:
        dfs: lista DataFrame'ów z danymi PM2.5.
    Returns:
        Połączony DataFrame z danymi PM2.5 oraz lista wspólnych stacji."""
    # Znajdź wspólne stacje
    return pd.concat(dfs, ignore_index=True, join="inner")

def read_data_from_csv(file_path: str) -> pd.DataFrame:
    """Wczytuje przetworzone dane z pliku CSV.
    Arguments:
        file_path: ścieżka do pliku CSV.
    Returns:
        DataFrame z danymi PM2.5."""
    df = pd.read_csv(file_path, header=[0,1])
    df.columns = pd.MultiIndex.from_tuples(
        [(a, "" if "Unnamed" in str(b) else b) for a, b in df.columns],
        names=["Miejscowość", "Kod stacji"]
    )
    df['Data'] = pd.to_datetime(df['Data'])
    return df

