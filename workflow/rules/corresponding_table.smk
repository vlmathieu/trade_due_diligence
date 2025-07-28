rule corresponding_table:
    output:
        'resources/inhouse/corresponding_table.csv'
    params:
        cmdCode         = config['uncomtrade']['cmdCode'],
        cmdDesc         = config['uncomtrade']['cmdDesc'],
        fao_items       = config['fao']['items'],
        fao_item_codes  = config['fao']['item_codes']
    log:
        'workflow/logs/corresponding_table.log'
    threads: 1
    conda:
        '../envs/polars.yaml'
    script:
        '../scripts/corresponding_table.py'