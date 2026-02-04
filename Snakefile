CONFIG_PATH = "config/task4.yaml"
configfile: CONFIG_PATH

YEARS = config["years"]

rule all:
    # Snakemake nie uruchomi rule'a jeśli pliki wyjściowe już istnieją
    input:
        # PM25
        expand("results/pm25/{year}/exceedance_days.csv", year=YEARS),
        expand("results/pm25/{year}/daily_means.csv", year=YEARS),
        expand("results/pm25/{year}/figures", year=YEARS),
        # PubMed
        expand("results/literature/{year}/pubmed_papers.csv", year=YEARS)

rule pm25_year:
    """
    Liczy statystyki PM2.5 dla jednego roku
    """
    input:
        config=CONFIG_PATH,
        load="src/pm25/load_data.py",
        analysis="src/pm25/data_analysis.py",
        viz="src/pm25/visualizations.py",
        runner="src/pm25/run_year.py"
    output:
        exceedance="results/pm25/{year}/exceedance_days.csv",
        daily="results/pm25/{year}/daily_means.csv",
        figs=directory("results/pm25/{year}/figures")
    params:
        year=lambda wc: int(wc.year),
        outdir=lambda wc: f"results/pm25/{wc.year}"
    log:
        "logs/pm25_{year}.log"
    shell:
        """
        mkdir -p {params.outdir}/figures
        python {input.runner} \
            --year {params.year} \
            --config {input.config} \
            > {log} 2>&1
        """

rule pubmed_year:
    """
    Pobiera dane literaturowe PubMed dla wybranego roku i słów kluczowych.
    """
    input:
        script="src/literature/pubmed_fetch.py",
        config=CONFIG_PATH
    output:
        csv="results/literature/{year}/pubmed_papers.csv"
    log:
        "logs/pubmed_{year}.log"
    shell:
        """
        python {input.script} \
            --year {wildcards.year} \
            --config {input.config} \
            > {log} 2>&1
        """