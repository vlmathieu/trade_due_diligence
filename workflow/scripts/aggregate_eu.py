from snakemake.script import snakemake
import logging
import polars as pl

# Log file edition
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

# Load data
input_data = pl.read_parquet(snakemake.input[0])
logging.info(f"\nInput data downloaded: \n {input_data.head(5)} \n")

# Transform EU country names and iso to "European Union" and "EU", respectively
eu = snakemake.params['eu']
logging.info(f"European countries considered: {eu}")
eu_iso = snakemake.params['eu_iso']
input_data_eu = (input_data
    .with_columns(
        reporter_desc = pl.when(pl.col('reporter_desc').is_in(eu))
                          .then(pl.lit('European Union'))
                          .otherwise(pl.col('reporter_desc'))
    )
    .with_columns(
        partner_desc = pl.when(pl.col('partner_desc').is_in(eu))
                         .then(pl.lit('European Union'))
                         .otherwise(pl.col('partner_desc'))
    )
    .with_columns(
        reporter_iso = pl.when(pl.col('reporter_iso').is_in(eu_iso))
                         .then(pl.lit('EU'))
                         .otherwise(pl.col('reporter_iso'))
    )
    .with_columns(
        partner_iso = pl.when(pl.col('partner_iso').is_in(eu_iso))
                        .then(pl.lit('EU'))
                        .otherwise(pl.col('partner_iso'))
    )
    .filter(pl.col('reporter_desc') != pl.col('partner_desc'))
)

# Sum netwgt and primary value by all other columns, i.e., 
# {'period', 'reporter_iso', 'reporter_desc', 'flow_code', 'partner_iso', 'partner_desc', 'cmd_code'} 
# ~ aggregation for duplicated trade flows EU <> partner | reporter <> EU
sum_cols = ['net_wgt', 'primary_value']

groupby_cols = [_ for _ in input_data_eu.columns if _ not in sum_cols]

input_data_eu = (
    input_data_eu
    .group_by(groupby_cols)
    .agg(pl.sum(sum_cols))
 )
logging.info(f"Data with european countries aggregated to EU:\n {input_data_eu.head(5)} \n {input_data_eu.describe()}")

# Save input data with aggregation of EU countries
input_data_eu.write_parquet(
    snakemake.output[0],
    compression='gzip'
    )
logging.info("Data saved.")