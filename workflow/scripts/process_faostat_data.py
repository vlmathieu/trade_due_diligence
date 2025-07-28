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

corr_table = pl.read_csv(
    source=snakemake.input[1],
    separator=';'
)
logging.info(f"\nCorresponding table downloaded: \n {corr_table} \n")

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
               .join(
                   corr_table.select(pl.exclude('faoDesc')),
                   left_on='Item Code',
                   right_on='faoCode')
)
logging.info(f"Filtered and joined data: \n {filter_data.head(5)} \n\n {filter_data.describe()} \n\n Columns: {filter_data.columns}\n")

# Process domestic consumption variables
cols = filter_data.select(pl.exclude(['Element', 'Value'])).columns
process_data = (filter_data
                .pivot(
                   on='Element',
                   index=cols,
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
logging.info(f"Proccessed data: \n {process_data.head(5)} \n\n {process_data.describe()} \n\n Columns: {process_data.columns}\n")

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