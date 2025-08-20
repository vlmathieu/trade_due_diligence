rule get_faostat:
    output:
        'resources/public/faostat.csv'
    params:
        fao_code = 'FO',
    log:
        'workflow/logs/get_faostat.log'
    threads: 1
    conda:
        '../envs/faostat.yaml'
    script:
        '../scripts/get_faostat.py'