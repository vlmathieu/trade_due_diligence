rule discrepancies_index:
    input:
        'results/intermediary/due_diligence_mirror_flows.csv'
    output:
        'results/intermediary/discrepancies_index.csv'
    log:
        'workflow/logs/discrepancies_index.log'
    threads: 4
    conda:
        '../envs/polars.yaml'
    script:
        '../scripts/discrepancies_index.py'