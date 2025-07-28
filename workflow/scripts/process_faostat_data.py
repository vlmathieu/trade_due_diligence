from snakemake.script import snakemake
import logging
import polars as pl

# Log file edition
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

# Load data
faostat_data = pl.read_csv(
    source=snakemake.input[0],
    separator=';',
    infer_schema_length=100000
)
logging.info(f"\nFAOSTAT data downloaded: \n {faostat_data.head(5)} \n")

# Filter data
logging.info(f"\nItem codes considered: {snakemake.params['item_codes']} \n")
logging.info(f"\nData cover trade from {snakemake.params['year_start']} to {snakemake.params['year_stop']-2}.\n")
filter_data = (faostat_data
               .filter(
                   
                   # Select item codes considered
                   pl.col('Item Code').is_in(snakemake.params['item_codes']),

                   # Keep data in specified time range
                   pl.col('Year') >= snakemake.params['year_start'],
                   pl.col('Year') <= snakemake.params['year_stop'] - 2,

                   # Remove non-country reporters
                   pl.col('Area Code (ISO3)').str.len_chars() == 3,
                   (~pl.col('Area Code (ISO3)')
                    .str.contains('|'.join(['XX', '_X', '\\d']))),

                   # Remove value in USD = keep volume data
                   ~pl.col('Unit').str.contains('USD')
                )
               .select(snakemake.params['col_keep'])
)

# Process domestic consumption variables
process_data = (filter_data
                .pivot(
                   on='Element',
                   index=['Area Code (ISO3)', 'Area', 'Item Code', 'Item', 'Year', 'Unit'],
                   values='Value'
               )
               .with_columns(
                   dom_consumption=(
                       pl.when(pl.col('Production') > 0)
                       .then(pl.col('Production')-pl.col('Export quantity'))
                       .otherwise(0)
                       ),
                   share_dom_consumption=(
                       pl.when(pl.col('Production') > 0)
                       .then((pl.col('Production')-pl.col('Export quantity'))/pl.col('Production')*100)
                       .otherwise(0)
                       )
               )
)
logging.info(f"Proccessed data: \n {process_data.head(5)} \n\n {process_data.describe()}")

# Save data
filter_data.write_csv(
        file=snakemake.output[0],
        separator=';'
        )
process_data.write_csv(
        file=snakemake.output[1],
        separator=';'
        )
logging.info("Data saved.")