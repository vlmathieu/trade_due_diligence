from snakemake.script import snakemake
import logging
import polars as pl

# Log file edition
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

# Build corresponding table
corr_table = pl.from_dict(
    {
        "cmdCode": snakemake.params['cmdCode'],
        "cmdDesc": snakemake.params['cmdDesc'],
        "faoCode": snakemake.params['fao_items'],
        "faoDesc": snakemake.params['fao_item_codes']
    }
)
logging.info(f"Corresponding table built:\n {corr_table}")

# Save corresponding table
corr_table.write_csv(
    file= snakemake.output[0],
    separator=';'
)
logging.info("Data saved.")
