DATA_DIR = "data/tables"
METADATA_FILE = "data/tables/metadata.csv"
PM25_RESULTS_DIR = "results/pm25"
CONFIG_PATH = "config/task4.yaml"
configfile: CONFIG_PATH

YEARS = config["years"]


rule all:
    """Sprawdza, czy pliki wyjściowe już istnieją"""
    input:
        # PM25
        expand(PM25_RESULTS_DIR + "/{year}/exceedance_days.csv", year=YEARS),
        expand(PM25_RESULTS_DIR + "/{year}/monthly_means_stations.csv", year=YEARS),
        expand(PM25_RESULTS_DIR + "/{year}/monthly_means_cities.csv", year=YEARS),
        expand(PM25_RESULTS_DIR + "/{year}/figures", year=YEARS),
        # PubMed
        expand("results/literature/{year}/pubmed_papers.csv", year=YEARS),
        # Raport
        "results/report_task4.md"

rule download_metadata:
    """Pobiera metadane."""
    input:
        script="src/pm25/download_metadata.py"
    output:
        csv=METADATA_FILE
    log:
        "logs/metadata.log"
    shell:
        """
        python {input.script} \
            --output {output.csv} \
            > {log} 2>&1
        """

rule download_data:
    """Pobiera dane dla wybranego roku."""
    input:
        meta=METADATA_FILE,
        script="src/pm25/download_data.py",
    output:
        csv = DATA_DIR + "/{year}.csv"
    params:
        year=lambda wc: wc.year,
        data_dir=DATA_DIR
    log:
        "logs/download_{year}.log"
    shell:
        """
        python {input.script} \
            --year {params.year} \
            --metadata {input.meta} \
            --output_dir {params.data_dir}
            > {log} 2>&1
        """

rule pm25_year:
    """
    Liczy statystyki PM2.5 dla jednego roku
    """
    input:
        data_csv = DATA_DIR + "/{year}.csv",
        metadata=METADATA_FILE,
        load="src/pm25/load_data.py",
        analysis="src/pm25/data_analysis.py",
        viz="src/pm25/visualizations.py",
        runner="src/pm25/run_year.py"
    output:
        exceedance = PM25_RESULTS_DIR + "/{year}/exceedance_days.csv",
        monthly_stations = PM25_RESULTS_DIR + "/{year}/monthly_means_stations.csv",
        monthly_cities = PM25_RESULTS_DIR + "/{year}/monthly_means_cities.csv",
        voivodeship = PM25_RESULTS_DIR + "/{year}/voivodeship_exceedance_days.csv",
        figs = directory(PM25_RESULTS_DIR + "/{year}/figures")
    params:
        year=lambda wc: int(wc.year),
        outdir=PM25_RESULTS_DIR,
        data_dir=DATA_DIR,
        config=CONFIG_PATH
    log:
        "logs/pm25_{year}.log"
    shell:
        """
        mkdir -p "{params.outdir}/{params.year}/figures"
        python {input.runner} \
            --year {params.year} \
            --config {params.config} \
            --metadata {input.metadata} \
            --data_dir {params.data_dir} \
            --output_dir {params.outdir}
            > {log} 2>&1
        """

rule pubmed_year:
    """
    Pobiera dane literaturowe PubMed dla wybranego roku i słów kluczowych.
    """
    input:
        script="src/literature/pubmed_fetch.py"
    output:
        csv="results/literature/{year}/pubmed_papers.csv"
    params:
        config=CONFIG_PATH
    log:
        "logs/pubmed_{year}.log"
    shell:
        """
        python {input.script} \
            --year {wildcards.year} \
            --config {params.config} \
            > {log} 2>&1
        """


rule report_task4:
    """
    Zbiera wyniki PM2.5 i PubMed dla wszystkich lat i generuje jeden wspólny raport.
    """
    input:
        pm25_files = expand(PM25_RESULTS_DIR + "/{year}/exceedance_days.csv", year=YEARS),
        pubmed_files = expand("results/literature/{year}/pubmed_papers.csv", year=YEARS),
        figures = expand(PM25_RESULTS_DIR + "/{year}/figures", year=YEARS),
        script = "src/report/generate_report.py",
        config = CONFIG_PATH
    output:
        report = "results/report_task4.md"
    params:
        years = YEARS,
        pm25_dir = PM25_RESULTS_DIR,
        lit_dir = "results/literature"
    log:
        "logs/report_task4.log"
    shell:
        """
        python {input.script} \
            --years {params.years} \
            --pm25_dir {params.pm25_dir} \
            --lit_dir {params.lit_dir} \
            --output {output.report} \
            --config {input.config} \
            > {log} 2>&1
        """