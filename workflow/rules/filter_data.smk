rule filter_data:
    input:
        'resources/public/uncomtrade_data.parquet.gzip'
    output:
        'results/input/input_data.parquet.gzip'
    params:
        year_stop       = config['years']['stop'],
        flow_keep       = config['flow_keep'],
        excluded_iso    = config['excluded_iso'],
        col_keep        = config['col_keep']
    log:
        'workflow/logs/filter_data.log'
    threads: 1
    conda:
        '../envs/polars.yaml'
    script:
        '../scripts/filter_data.py'