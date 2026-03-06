rule get_uncomtrade:
    input:
        'resources/inhouse/due_diligence_codes_hs6.csv',
        'resources/inhouse/placebo_codes_hs6.csv'
        # 'resources/inhouse/products_under_due_diligence.json'
    output:
        'resources/public/due_diligence_data.parquet.gzip',
        'resources/public/placebo_data.parquet.gzip'
    params:
        year_start      = config['years']['start'],
        year_stop       = config['years']['stop'],
        flowCode        = list(str(flow) for flow in config['flowCode']),
        apikey          = os.environ['comtrade_apikey']
    log:
        'workflow/logs/get_uncomtrade.log'
    threads: 1
    conda:
        '../envs/comtradeapicall.yaml'
    script:
        '../scripts/get_uncomtrade.py'