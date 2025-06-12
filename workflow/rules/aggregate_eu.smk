rule aggregate_eu:
    input:
        'results/input/input_data.parquet.gzip'
    output:
        'results/intermediary/input_data_eu.parquet.gzip'
    params:
        eu      = config['eu']['countries'],
        eu_iso  = config['eu']['iso']
    log:
        'workflow/logs/aggregate_eu.log'
    threads: 1
    conda:
        '../envs/polars.yaml'
    script:
        '../scripts/aggregate_eu.py'