# Analiza jakości powietrza (PM2.5) i trendów publikacji

Projekt służy do automatycznej analizy danych o zanieczyszczeniu powietrza (PM2.5) w Polsce oraz powiązanych trendów w publikacjach naukowych (PubMed). Całość procesu jest zarządzana przez Snakemake, co gwarantuje powtarzalność i optymalizację obliczeń.

## Struktura projektu

```text
├── config/
│   └── task4.yaml              # Konfiguracja (lata, słowa kluczowe, e-mail)
├── data/
│   ├── gios_ids.json           # Mapowanie ID stacji GIOS
│   └── wojewodztwa-min.geojson # Dane geograficzne do map
├── src/
│   ├── literature/             # Skrypty do pobierania danych z PubMed
│   ├── pm25/                   # Skrypty do pobierania i analizy danych PM2.5
│   └── report/                 # Skrypt generujący końcowy raport Markdown
├── results/                    # Katalog z wynikami (generowany automatycznie)
├── tests/                      # Testy jednostkowe (pytest)
├── Snakefile                   # Definicja reguł pipeline'u
├── requirements.txt            # Zależności Python
└── README.md

```

## Konfiguracja (`config/task4.yaml`)

Działanie pipeline'u jest sterowane przez plik YAML. Pozwala on na łatwą zmianę zakresu analizy bez ingerencji w kod.

**Przykładowa konfiguracja:**

```yaml
years: [2022]                                   # Lista lat do analizy (PM2.5 i PubMed)
cities: ["Warszawa", "Katowice"]                # Miasta, dla których generowane są szczegółowe wykresy trendów
entrez_email: "test@example.com"                # Adres e-mail wymagany przez API NCBI (PubMed)
max_records: 100                                # Maksymalna liczba pobieranych artykułów na rok
journals: ["Nature", "Science"]                 # Lista czasopism do filtrowania zapytań
keywords: ["Air Pollution", "PM2.5"]            # Słowa kluczowe wyszukiwania

```

## Instalacja

Zainstaluj wymagane biblioteki:
```bash
pip install -r requirements.txt

```

## Uruchamianie i przykładowy scenariusz działania

Pipeline jest sterowany przez plik konfiguracyjny `config/task4.yaml`. Poniżej przedstawiony został przykładowy scenariusz uruchamiania pipeline'u.

### Krok 1: Pierwsze uruchomienie

Użytkownik ustawia w `config/task4.yaml`:

```yaml
years: [2021, 2024]

```

Następnie uruchamia pipeline:

```bash
snakemake --cores 1

```

**Rezultat:**

* System pobiera i przetwarza dane PM2.5 dla lat 2021 i 2024.
* System pobiera dane z PubMed dla lat 2021 i 2024.
* Generowany jest raport `results/report_task4.md` obejmujący lata {2021, 2024}.

### Krok 2: Zmiana konfiguracji 

Użytkownik zmienia w `config/task4.yaml`:

```yaml
years: [2019, 2024]

```

Uruchamia ponownie tę samą komendę:

```bash
snakemake --cores 1

```

**Oczekiwane zachowanie:**
Pipeline wykonuje **tylko brakujące kroki**:

1. Liczy PM2.5 oraz pobiera dane PubMed **tylko dla roku 2019**.
2. **Pomija** ponowne liczenie roku 2024 (artefakty już istnieją).
3. Generuje nowy raport zbiorczy `results/report_task4.md` dla lat {2019, 2024}.

### Weryfikacja

Poprawność działania (brak ponownego przeliczania roku 2024) jest weryfikowana poprzez analizę logów Snakemake: Dla reguł dotyczących roku 2024 (np. `pm25_year`, `pubmed_year`) nie pojawia się status "Running", a jedynie dla roku 2019. Snakemake zgłosi wykonanie zadań tylko dla nowych danych oraz reguły `report_task4` (ponieważ zmieniły się wejścia).

### Zaawansowane opcje uruchamiania

Domyślnie Snakemake decyduje o ponownym uruchomieniu zadań na podstawie kilku reguł ([link do dokumnetacji](https://snakemake.readthedocs.io/en/v7.14.1/executing/cli.html)): 

* `mtime`: czas modyfikacji plików wejściowych (czy którykolwiek z nich ma nowszy mtime niż plik wyjściowy);
* `params`: zmiana parametrów z sekcji *params*;
* `input`: zmiana listy logicznych wejść reguły;
* `software-env`: zmiana w środowisku wykonania;
* `code`: zmiany w kodzie reguły.

Zapewnia to spójność wyników z aktualnym kodem workflow, konfiguracją oraz środowiskiem uruchomieniowym. Można modyfikować triggery flagą `--rerun-triggers`. Przykładowo, uruchomienie pipeline'u komendą

```bash
snakemake --cores 1 --rerun-triggers mtime
```

W takim trybie Snakemake ponownie uruchamia zadania tylko na podstawie czasu modyfikacji plików.


## Wyniki

Wyniki są zapisywane w katalogu `results/` z podziałem na lata, co zapobiega nadpisywaniu danych:

* `results/pm25/{ROK}/` – dane CSV, wykresy i mapy dla danego roku.
* `results/literature/{ROK}/` – dane CSV z publikacjami.
* `results/report_task4.md` – zbiorczy raport w formacie Markdown.

## Testowanie

W projekcie zaimplementowano testy jednostkowe przy użyciu `pytest`, które pokrywają kluczowe funkcjonalności przetwarzania i analizy danych. Sprawdzają one m.in. poprawność parsowania dat z PubMed, mapowanie kodów stacji oraz analizy statystyczne.

Uruchomienie testów:

```bash
pytest

```
