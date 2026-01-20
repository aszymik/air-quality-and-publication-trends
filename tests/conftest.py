import pandas as pd
import pytest

@pytest.fixture
def sample_df():
    data = {
        ('Data',''): pd.to_datetime([
            '2022-01-01', '2022-01-02',
            '2022-02-01', '2022-02-02'
        ]),
        ('Warszawa','stacja1'): [10, None, 30, 40],
        ('Warszawa','stacja2'): [20, 30, 40, 50],
        ('Kraków','stacja1'): [5, 15, 25, None],
    }
    df = pd.DataFrame(data)
    df.columns = pd.MultiIndex.from_tuples(
        df.columns,
        names=["Miejscowość", "Kod stacji"]
    )
    return df

@pytest.fixture
def midnight():
    return pd.DataFrame({
        'Data': ['2021-01-15 14:30:00', '2021-01-15 00:00:00', '2023-01-17 23:59:00'],
        'Temperature': [5, -2, 0]
    })

@pytest.fixture
def sample_metadata():
    return pd.DataFrame({
        'Kod stacji': ['stacja1', 'stacja2', 'stacja3'],
        'Stary kod stacji': ['old1, old2', pd.NA, 'old3'],
        'Miejscowość': ['Warszawa', 'Kraków', 'Gdańsk'],
        'Województwo': ['Mazowieckie', 'Małopolskie', 'Pomorskie'],
    })

@pytest.fixture
def df_without_multiindex():
    data = {
        'Data': pd.to_datetime([
            '2022-01-01', '2022-01-02',
            '2022-02-01', '2022-02-02'
        ]),
        'stacja1': [10, None, 30, 40],
        'stacja2': [20, 30, 40, 50],
        'stacja3': [5, 15, 25, None],
    }
    return pd.DataFrame(data)

@pytest.fixture
def code_to_city_dict():
    return {
        'stacja1': 'Warszawa',
        'stacja2': 'Kraków',
    }

@pytest.fixture
def code_to_voivodeship_dict():
    return {
        'stacja1': 'Mazowieckie',
        'stacja2': 'Małopolskie',
    }
