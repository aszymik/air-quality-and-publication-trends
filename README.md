# Analiza danych dotyczących poziomu zanieczyszczenia PM2.5


Na podstawie danych udostępnionych przez GIOŚ o jakości powietrza w Polsce projekt przeprowadza analizę wraz z wykresami, istotnymi porównaniami oraz wykrywaniem dni z przekroczeniem norm na przestrzeni lat (>15 µg/m³). 

Skupia się na pomiarach stężeń drobnego pyłu PM2.5, ze względu na wyjątkową szkodliwość wdychania tego pyłu na układ oddechowy i krwioobieg.

Źródło danych: https://powietrze.gios.gov.pl/

Źródło mapy używanej do obsługi granic województw: https://github.com/ppatrzyk/polska-geojson

## Wymagania
Python 3.11

Wymagane biblioteki zostały wymienione w pliku **requirements.txt**
Główne:
- pandas (przetrwarzanie danych)
- matplotlib (wizualizacja)
- requests (pobieranie dnaych z API GIOŚ)
- pytest (testy)

## Instalacja
```bash
git clone https://github.com/aszymik/air-quality-and-publication-trends.git

cd polish-air-quality-trends

pip install -r requirements.txt
```

## Użycie

Najlepiej uruchomić gotową analizę w notatniku:

```bash

jupyter notebook main.ipynb
```
Notatnik zawiera pełną wizualizację analizy.

### Przykład użycia w kodzie
```python
from scripts.load_data import read_data_from_csv
from scripts.data_analysis import  get_chosen_monthly_means, get_who_norm_exceeding_days

df = read_data_from_csv('data/pm25_gios_2015_2018_2021_2024.csv')

monthly_chosen = get_chosen_monthly_means(df, chosen_years=[2015, 2024], chosen_cities=['Warszawa', 'Katowice'])

exceeding_days = get_who_norm_exceeding_days(df)
```
## Struktura projektu
    polish-air-qaulity-trends/

        data/
            wojewodztwa-min.geojson

        scripts/
            data_analysis.py
            load_data.py
            visualizations.py
            
        tests/
            conftest.py
            test_data_analysis.py
            test_load_data.py

        main.ipynb
        README.md
        requirements.txt

## Uruchamianie testów
```bash
pytest

```
## Autorzy 
Maja Kończak

Anna Szymik


