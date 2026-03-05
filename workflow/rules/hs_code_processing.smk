rule hs_code_processing:
    input:
        'resources/public/hs_classification.csv',
        'resources/inhouse/due_diligence_codes.csv',
        'resources/inhouse/placebo_codes.csv'
    output:
        'resources/inhouse/due_diligence_codes_hs6.csv',
        'resources/inhouse/placebo_codes_hs6.csv'
    log:
        'workflow/logs/hs_code_processing.log'
    threads: 4
    conda:
        '../envs/polars.yaml'
    script:
        '../scripts/hs_code_processing.py'