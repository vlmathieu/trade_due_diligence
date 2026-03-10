rule mirror_flow:
    input:
        'resources/public/due_diligence_data.parquet.gzip'
    output:
        'results/intermediary/due_diligence_mirror_flows.csv'
    params:
        year_stop       = config['years']['stop'],
        excluded_iso    = config['excluded_iso'],
        col_keep        = config['col_keep']['uncomtrade']
    log:
        'workflow/logs/mirror_flow.log'
    threads: 4
    conda:
        '../envs/polars.yaml'
    script:
        '../scripts/mirror_flow.py'