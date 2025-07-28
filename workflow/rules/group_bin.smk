rule group_bin:
    input:
        'results/input/input_uncomtrade.parquet.gzip',
    output:
        'results/intermediary/uncomtrade_group_bin.parquet.gzip'
    params:
        policies   = config['policy'],
        risky      = config['risky']['countries']
    log:
        'workflow/logs/group_bin.log'
    threads: 1
    conda:
        '../envs/polars.yaml'
    script:
        '../scripts/group_bin.py'