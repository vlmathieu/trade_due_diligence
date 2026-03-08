rule filter_data:
    input:
        'resources/public/due_diligence_data.parquet.gzip',
        'resources/public/placebo_data.parquet.gzip'
    output:
        'results/input/due_diligence_input.csv',
        'results/input/placebo_input.csv'
    params:
        year_stop       = config['years']['stop'],
        flow_keep       = config['flow_keep'],
        excluded_iso    = config['excluded_iso'],
        col_keep        = config['col_keep']['uncomtrade']
    log:
        'workflow/logs/filter_data.log'
    threads: 1
    conda:
        '../envs/polars.yaml'
    script:
        '../scripts/filter_data.py'