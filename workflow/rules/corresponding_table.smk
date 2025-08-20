rule corresponding_table:
    output:
        'resources/inhouse/corresponding_table.csv'
    log:
        'workflow/logs/corresponding_table.log'
    threads: 1
    conda:
        '../envs/polars.yaml'
    script:
        '../scripts/corresponding_table.py'