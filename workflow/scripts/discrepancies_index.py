from snakemake.script import snakemake
import logging
import polars as pl

# Log file edition
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
    
# Load mirror flow data
mirror_flows = pl.read_csv(snakemake.input[0], separator=',')

# Compute discrepancies index
discrepancies_index = (
    mirror_flows
    .with_columns(pl.col('cmd_code').cast(pl.String))
    .filter(pl.col('cmd_code').str.starts_with('4403'))
    .with_columns(
        pl.col('primary_value_exp').fill_null(0),
        pl.col('primary_value_imp').fill_null(0)
    )
    .with_columns(
        discrepancies_index=(
            (pl.col('primary_value_exp') - pl.col('primary_value_imp')) /
            (pl.col('primary_value_exp') + pl.col('primary_value_imp'))
            )
    )
    .select(['exporter_iso', 'discrepancies_index'])
    .group_by('exporter_iso', maintain_order=True)
    .mean()
    .sort('discrepancies_index', descending=True)
)
logging.info(f"\nDiscrepancies index:\n {discrepancies_index}\n")

# Save discrepancies index
discrepancies_index.write_csv(
    snakemake.output[0],
    separator=','
    )
