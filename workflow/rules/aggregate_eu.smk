rule aggregate_eu:
    input:
        'results/input/input_uncomtrade.parquet.gzip'
    output:
        'results/input/input_uncomtrade_eu.parquet.gzip'
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