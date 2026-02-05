import argparse
import pandas as pd
import os
import sys
import re
import matplotlib.pyplot as plt
from collections import Counter
from pathlib import Path


STOPWORDS = set([
    "the", "of", "and", "in", "to", "a", "for", "on", "with", "is", "by", "an",
    "at", "from", "as", "validation", "using", "study", "analysis", "during",
    "between", "impact", "effect", "effects", "associated", "association",
    "concentration", "concentrations", "levels", "quality", "pollution", "particulate", "matter"
])

def load_pm25_data(years, pm25_dir):
    dfs = []
    print(f"DEBUG: Szukam danych PM2.5 w {pm25_dir} dla lat {years}")
    
    for year in years:
        file_path = Path(pm25_dir) / str(year) / "exceedance_days.csv"
        if file_path.exists():
            try:
                df = pd.read_csv(file_path, index_col=0)
                if df.empty: continue
                df.columns = df.columns.astype(str)
                dfs.append(df)
            except Exception as e:
                print(f"Błąd przy wczytywaniu {file_path}: {e}")
    
    if not dfs: return pd.DataFrame()

    try:
        return pd.concat(dfs, axis=1, join='outer').sort_index()
    except Exception as e:
        print(f"Błąd łączenia: {e}")
        return pd.DataFrame()

def load_literature_data(years, lit_dir):
    all_papers = []
    for year in years:
        file_path = Path(lit_dir) / str(year) / "pubmed_papers.csv"
        if file_path.exists():
            try:
                df = pd.read_csv(file_path)
                if not df.empty: all_papers.append(df)
            except: pass
    
    if not all_papers: return pd.DataFrame()
    return pd.concat(all_papers, ignore_index=True)

def analyze_titles_and_plot(lit_df, years, output_dir):
    """
    Tworzy wykres porównawczy najczęstszych słów dla pierwszego i ostatniego roku.
    Zapisuje plik PNG i zwraca ścieżkę względną do wstawienia w MD.
    """
    if lit_df.empty or 'title' not in lit_df.columns:
        return None

    # Wybieramy pierwszy i ostatni rok do porównania
    year_start = min(years)
    year_end = max(years)
    
    compare_years = [year_start]
    if year_end != year_start:
        compare_years.append(year_end)

    # Zbieranie słów per rok
    word_counts = {}
    
    for y in compare_years:
        titles = lit_df[lit_df['year'] == y]['title'].dropna().astype(str).tolist()
        words = []
        for t in titles:
            # Tokenizacja: małe litery, tylko znaki alfanumeryczne
            tokens = re.findall(r'\b[a-z]{3,}\b', t.lower())
            # Filtrowanie po stopwords
            tokens = [w for w in tokens if w not in STOPWORDS]
            words.extend(tokens)
        word_counts[y] = Counter(words)

    if not any(word_counts.values()):
        return None

    # Top 10 słów dla każdego z wybranych lat
    all_top_words = set()
    for y in compare_years:
        for word, _ in word_counts[y].most_common(10):
            all_top_words.add(word)
    
    sorted_words = sorted(list(all_top_words))
    
    # Dataframe do wykresu
    plot_data = {'Słowo': sorted_words}
    for y in compare_years:
        plot_data[str(y)] = [word_counts[y].get(w, 0) for w in sorted_words]
    
    df_plot = pd.DataFrame(plot_data).set_index('Słowo')

    plt.figure(figsize=(10, 6))
    df_plot.plot(kind='bar', width=0.8, color=['skyblue', 'salmon'][:len(compare_years)])
    plt.title(f"Najczęstsze słowa w tytułach: {year_start} vs {year_end}")
    plt.ylabel("Liczba wystąpień")
    plt.xlabel("Słowo")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    filename = "title_analysis_comparison.png"
    output_path = Path(output_dir) / filename
    plt.savefig(output_path)
    plt.close()
    
    return filename

def add_visualizations(md, years, pm25_dir):
    md.append("## Wizualizacje PM2.5\n")
    md.append("Poniżej przedstawiono wybrane wykresy wygenerowane w procesie analizy.\n")

    for year in years:
        md.append(f"### Rok {year}")
        figs_path = Path(pm25_dir) / str(year) / "figures"
        images_to_show = [
            ("voivodeship_exceeding_days_map.png", "Mapa przekroczeń w województwach"),
            ("who_exceeding_days.png", "Dni z przekroczeniem normy WHO (3 najlepsze i 3 najgorsze stacje)"),
            ("trends_chosen_cities.png", "Trendy dla wybranych miast")
        ]
        
        found_any = False
        for filename, description in images_to_show:
            img_full_path = figs_path / filename
            
            if img_full_path.exists():
                # Markdown potrzebuje względnej ścieżki
                rel_path = str(img_full_path).replace("results/", "").replace("results\\", "")
                md.append(f"**{description}**")
                md.append(f"![{description}]({rel_path})\n")
                found_any = True
        
        if not found_any:
            md.append("_Brak wygenerowanych wykresów dla tego roku._\n")
        
        md.append("---\n")
    return md

def generate_markdown(years, pm25_df, lit_df, pm25_dir, output_file_path):
    md = []
    md.append(f"# Raport Task 4: Jakość powietrza i analiza literatury\n")
    md.append(f"**Zakres lat:** {', '.join(map(str, years))}\n")

    # --- PM2.5 ---
    md.append("## 1. Analiza PM2.5: Dni z przekroczeniem normy WHO")
    if not pm25_df.empty:
        try:
            pm25_numeric = pm25_df.apply(pd.to_numeric, errors='coerce').fillna(0)
            pm25_df['Suma'] = pm25_numeric.sum(axis=1)
            top_stations = pm25_df.sort_values('Suma', ascending=False).head(10)
            display_df = top_stations.drop(columns=['Suma'])
            
            md.append("### Top 10 stacji z największą liczbą przekroczeń (wszystkie lata)")
            try:
                md.append(display_df.fillna("-").to_markdown())
            except ImportError:
                md.append(display_df.fillna("-").to_string())
            md.append("\n*Legenda: '-' oznacza brak pomiarów.*")
        except Exception as e:
            md.append(f"Błąd tabeli: {e}")
    else:
        md.append("_Brak danych PM2.5._")

    md.append("\n---\n")

    # Wizualizacje
    md = add_visualizations(md, years, pm25_dir)
    md.append("\n---\n")

    # --- Literatura ---
    md.append("## 2. Analiza Literatury (PubMed)")
    if not lit_df.empty:
        md.append(f"**Łączna liczba publikacji:** {len(lit_df)}\n")
        
        # Analiza słów kluczowych
        md.append("### Porównanie tematów (najczęstsze słowa w tytułach)")
        output_dir = os.path.dirname(output_file_path)
        plot_filename = analyze_titles_and_plot(lit_df, years, output_dir)
        
        if plot_filename:
            md.append(f"Wykres przedstawia najczęstsze słowa (z wyłączeniem stopwords) w tytułach publikacji dla lat: {min(years)} i {max(years)}.\n")
            md.append(f"![Analiza słów]({plot_filename})\n")
        else:
            md.append("_Brak danych do analizy tekstowej (zbyt mało tytułów)._\n")

        # Trend liczby publikacji
        try:
            md.append("### Trend liczby publikacji")
            trend = lit_df['year'].value_counts().sort_index().reset_index()
            trend.columns = ['Rok', 'Liczba']
            try:
                md.append(trend.to_markdown(index=False))
            except ImportError:
                md.append(trend.to_string(index=False))
            
            md.append("\n### Przykładowe publikacje")
            sample = lit_df.sample(min(5, len(lit_df)))
            for _, row in sample.iterrows():
                title = str(row.get('title', 'Brak tytułu')).replace('\n', ' ')
                md.append(f"- {title}")
        except Exception as e:
            md.append(f"Błąd sekcji literatury: {e}")
    else:
        md.append("_Brak publikacji._")

    md.append("\n---\n")

    return "\n".join(md)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--years", nargs="+", type=int, required=True)
    parser.add_argument("--pm25_dir", required=True)
    parser.add_argument("--lit_dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    if not os.path.exists(args.pm25_dir):
        print(f"Błąd: katalog {args.pm25_dir} nie istnieje!")
        sys.exit(1)

    pm25_data = load_pm25_data(args.years, args.pm25_dir)
    lit_data = load_literature_data(args.years, args.lit_dir)
    report_content = generate_markdown(args.years, pm25_data, lit_data, args.pm25_dir, args.output)

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as f:
        f.write(report_content)
    print(f"Raport zapisany: {args.output}")

if __name__ == "__main__":
    main()