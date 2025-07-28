rule get_uncomtrade:
    output:
        'resources/public/uncomtrade.parquet.gzip'
    params:
        year_start      = config['years']['start'],
        year_stop       = config['years']['stop'],
        cmdCode         = config['uncomtrade']['cmdCode'],
        flowCode        = list(str(flow) for flow in config['flowCode']),
        apikey          = os.environ['comtrade_apikey']
    log:
        'workflow/logs/get_uncomtrade.log'
    threads: 1
    conda:
        '../envs/comtradeapicall.yaml'
    script:
        '../scripts/get_uncomtrade_data.py'