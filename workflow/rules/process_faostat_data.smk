rule process_faostat_data:
    input:
        'resources/public/faostat.csv',
        'resources/inhouse/corresponding_table.csv'
    output:
        'results/input/input_faostat.csv',
        'results/intermediary/processed_faostat.csv'
    params:
        year_start      = config['years']['start'],
        year_stop       = config['years']['stop'],
        excluded_iso    = config['excluded_iso'],
        col_keep        = config['col_keep']['fao']
    log:
        'workflow/logs/process_faostat_data.log'
    threads: 1
    conda:
        '../envs/polars.yaml'
    script:
        '../scripts/process_faostat_data.py'