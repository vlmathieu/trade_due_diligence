from snakemake.script import snakemake
import logging
import polars as pl
import faostat
import time

# Log file edition
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

# Get parameter
fao_code = snakemake.params['fao_code']

# Collect parameters from faostat database
## Elements: production, export, import
elements = pl.from_dataframe(faostat.get_par_df(fao_code, 'element'))
elements_code = (elements
                .select(pl.col('code'))
                .to_series()
                .to_list()
)

## Items: wood products listed in database
items = pl.from_dataframe(faostat.get_par_df(fao_code, 'item'))
items_code = (items
            .filter(~pl.col('code').str.contains('>'))
            .select('code')
            .to_series()
            .to_list()
)

# Collect data
logging.info("Starting to download...\n")
start_time = time.time()
faostat_data = pl.concat([
    pl.from_dataframe(
        faostat.get_data_df(
            code    = fao_code,
            pars    = {'element':elements_code, 'item':item},
            coding  = {'area': 'ISO3'}
        )
    )
    for item in items_code
])
logging.info(f"Downloading data took {round((time.time() - start_time), 2)} seconds\n")

logging.info(f"\nDataframe head: \n {faostat_data.head(5)} \n")
logging.info(f"\nDataframe size (rows, columns): {faostat_data.shape} \n")

# Save data if check list passed
check_list = (
    sorted(items_code) == sorted(faostat_data.select(pl.col("Item Code")).unique().to_series().to_list())
)
if check_list:
    logging.info('Data have been checked.\n')   
    faostat_data.write_csv(
        file=snakemake.output[0],
        separator=';'
        )
else:
    logging.info('Issues found in data download.\n')
