rule products_under_due_diligence:
    output:
        'resources/inhouse/products_under_due_diligence.json'
    log:
        'workflow/logs/products_under_due_diligence.log'
    threads: 1
    conda:
        '../envs/polars.yaml'
    script:
        '../scripts/products_under_due_diligence.py'