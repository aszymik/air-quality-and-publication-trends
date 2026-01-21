configfile: "config/task4.yaml"

YEARS = config["years"]

rule all:
    # Snakemake nie uruchomi rule'a jeśli pliki output już istnieją
    input:
        expand("results/pm25/{year}/exceedance_days.csv", year=YEARS),
        expand("results/pm25/{year}/daily_means.csv", year=YEARS)

rule pm25_year:
    """
    Liczy statystyki PM2.5 dla jednego roku
    """
    input:
        config="config/task4.yaml",
        load="scripts/pm25/load_data.py",
        analysis="scripts/pm25/data_analysis.py",
        viz="scripts/pm25/visualizations.py",
        runner="scripts/pm25/run_year.py"
    output:
        exceedance="results/pm25/{year}/exceedance_days.csv",
        daily="results/pm25/{year}/daily_means.csv"
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
