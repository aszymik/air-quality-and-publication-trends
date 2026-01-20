import pandas as pd
from scripts.load_data import (
    get_code_mappings,
    rename_columns,
    add_multiindex,
    change_midnight_measurements
    )

def test_rename_columns():
    df = pd.DataFrame({
        'Data': [1, 2],
        'old1': [3, 4],
        'old2': [5, 6],
        'new3': [7, 8]
    })
    old_to_new_code = {
        'old1': 'A',
        'old2': 'B'
    }
    result = rename_columns(df, old_to_new_code)
    assert list(result.columns) == ['Data', 'A', 'B', 'new3']

def test_change_midnight_measurements(midnight):
    result = change_midnight_measurements(midnight)
 
    assert result.loc[0, 'Data'] == pd.Timestamp('2021-01-15 14:30:00')
    assert result.loc[1, 'Data'] == pd.Timestamp('2021-01-14 23:59:59')
    assert result.loc[2, 'Data'] == pd.Timestamp('2023-01-17 23:59:00')

def test_change_midnight_measurements_none(midnight):
    midnight['Data'] = [pd.NaT, '2023-01-15 14:30:00', '2023-01-17 23:59:59']
    result = change_midnight_measurements(midnight)
    print(type(result.loc[0, 'Data']))
    # assert pd.isna(result.loc[0, 'Data'])
    assert result.loc[1, 'Data'] == pd.Timestamp('2023-01-15 14:30:00')

def test_get_code_mappings(sample_metadata):
    old_to_new, code_to_city, code_to_voivodeship = get_code_mappings(sample_metadata)

    assert old_to_new == {
        'old1': 'stacja1',
        'old2': 'stacja1',
        'old3' : 'stacja3'
    }
    assert code_to_city == {
        'stacja1': 'Warszawa',
        'stacja2': 'Kraków',
        'stacja3': 'Gdańsk'
    }
    assert code_to_voivodeship == {
        'stacja1': 'Mazowieckie',
        'stacja2': 'Małopolskie',
        'stacja3': 'Pomorskie'
    }

def test_get_code_mappings_whitespaces():
    metadata = pd.DataFrame({
        'Kod stacji': ['stacja1'],
        'Stary kod stacji': [' old1, old2 '],
        'Miejscowość': [' Warszawa '],
        'Województwo': [' Mazowieckie ']
    })
    old_to_new, _, code_to_voivodeship = get_code_mappings(metadata)
    assert 'old1' in old_to_new
    assert 'old2' in old_to_new
    assert old_to_new['old1'] == 'stacja1'
    assert ' old1 ' not in old_to_new
    assert code_to_voivodeship['stacja1'].strip() == 'Mazowieckie'

def test_get_code_mappings_no_old_codes(sample_metadata):

    sample_metadata['Stary kod stacji'] = [pd.NA, pd.NA, pd.NA]

    old_to_new, code_to_city, _ = get_code_mappings(sample_metadata)
    assert len(old_to_new) == 0
    assert code_to_city['stacja1'] == 'Warszawa'
    assert code_to_city['stacja2'] == 'Kraków'

def test_add_multiindex(df_without_multiindex, code_to_city_dict):
    result = add_multiindex(df_without_multiindex, code_to_city_dict)
    assert ('Warszawa','stacja1') in result.columns
    assert ('Kraków','stacja2') in result.columns
    assert ('Nieznane','stacja3') in result.columns
