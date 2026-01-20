import pandas as pd
from scripts.data_analysis import (
    get_monthly_means_for_stations,
    get_chosen_monthly_means,
    get_monthly_means_for_cities,
    get_who_norm_exceeding_days,
    get_max_and_min_k_stations,
    get_voivodeship_exceeding_days,
)

def test_get_monthly_means_for_stations(sample_df):
    result = get_monthly_means_for_stations(sample_df)
    assert len(result) == 2

    jan_mean = 10.0  # styczeń 2022
    assert result.loc[(2022, 1), ('Warszawa','stacja1')] == jan_mean
    assert result.loc[(2022, 1), ('Warszawa','stacja2')] == 25.0
    assert result.loc[(2022, 1), ('Kraków','stacja1')] == 10.0  
    assert result.loc[(2022, 2), ('Warszawa','stacja1')] == 35.0

def test_get_monthly_means_for_stations_none(sample_df):
    sample_df[('Warszawa','stacja1')] = [None, None, None, None]
    sample_df[('Kraków','stacja1')] = [None, None, None, None]
    result = get_monthly_means_for_stations(sample_df)

    assert pd.isna(result.loc[(2022, 1), ('Warszawa','stacja1')])
    assert pd.isna(result.loc[(2022, 1), ('Kraków','stacja1')])

def test_get_chosen_monthly_means(sample_df):
    result = get_chosen_monthly_means(
        sample_df,
        chosen_years=[2022],
        chosen_cities=['Warszawa']
    )

    assert len(result) == 2  # 2 miesiące
    assert set(result['Miejscowość']) == {'Warszawa'}

    jan_value = result.loc[result['Miesiąc'] == 1, 'PM2.5'].iloc[0]  # średnia za styczeń
    assert jan_value == (10 + 20 + 30) / 3

def test_get_monthly_means_for_cities(sample_df):
    result = get_monthly_means_for_cities(sample_df)

    assert len(result) == 2
    assert 'Rok' in result.columns
    assert 'Miesiąc' in result.columns
    assert result.loc[(result['Rok'] == 2022) & (result['Miesiąc'] == 1), 'Warszawa'].iloc[0] == ((10+20)/2 + 30) / 2
    assert result.loc[(result['Rok'] == 2022) & (result['Miesiąc'] == 2), 'Warszawa'].iloc[0] == ((30+40)/2 + (50+40)/2) / 2

def test_get_who_norm_exceeding_days(sample_df):
    result = get_who_norm_exceeding_days(sample_df)

    assert 2022 in result.columns
    assert result.loc[('Warszawa','stacja2'), 2022] == 4

def test_get_max_and_min_k_stations():

    data = {
        2022: [5, 10, 2, 20]
    }
    df = pd.DataFrame(
        data,
        index=['A', 'B', 'C', 'D']
    )
    result = get_max_and_min_k_stations(df, chosen_year=2022, k=1)
    assert list(result.index) == ['C', 'D']

def test_get_voivodeship_exceeding_days():
    df = pd.DataFrame({
        ('Data',''): pd.to_datetime(['2022-01-01', '2022-01-02', '2022-01-03']),
        ('Warszawa','stacja1'): [10, 20, 5],
        ('Warszawa','stacja2'): [None, 5, 30],
        ('Kraków','stacja3'): [16, 14, 10],
    })
    df.columns = pd.MultiIndex.from_tuples(
        df.columns,
        names=["Miejscowość", "Kod stacji"]
    )
    code_to_voivodeship = {
        'stacja1': 'Mazowieckie',
        'stacja2': 'Mazowieckie',
        'stacja3': 'Małopolskie',
    }
    result = get_voivodeship_exceeding_days(df, code_to_voivodeship, threshold=15)

    assert 2022 in result.columns
    assert result.loc['Mazowieckie', 2022] == 2
    assert result.loc['Małopolskie', 2022] == 1
