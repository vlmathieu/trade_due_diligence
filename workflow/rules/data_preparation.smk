rule data_preparation:
    input:
        'resources/public/due_diligence_data.parquet.gzip',
        'resources/public/placebo_data.parquet.gzip',
        'resources/inhouse/due_diligence_codes_hs6.csv',
        'resources/inhouse/placebo_codes_hs6.csv'
    output:
        'results/input/due_diligence_input.csv',
        'results/input/placebo_input.csv'
    params:
        year_stop       = config['years']['stop'],
        flow_keep       = config['flow_keep'],
        excluded_iso    = config['excluded_iso'],
        col_keep        = config['col_keep']['uncomtrade']
    log:
        'workflow/logs/data_preparation.log'
    threads: 1
    conda:
        '../envs/polars.yaml'
    script:
        '../scripts/data_preparation.py'