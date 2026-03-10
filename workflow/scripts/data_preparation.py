from snakemake.script import snakemake
import logging
import polars as pl
import re

def format_col_names(col_names: list) -> dict:
    '''
    Function that returns a mapping between column names and their formated 
    version as a dictionnary.

    Parameters
    ----------
    col_names : list of strings
        The list of column names to format.

    Returns
    -------
    col_mapping : dictionnary
        Mapping between column names and their formated version.
    '''

    # Insert an underscore before a uppercase
    col_names_formated = [re.sub(r"([a-z])([A-Z])", "\\1_\\2", s) for s in col_names]

    # Replace space by underscore
    col_names_formated = [re.sub(" ", "_", s) for s in col_names_formated]

    # Put every column name to lowercase
    col_names_formated = [s.lower() for s in col_names_formated]

    # Create mapping between column names and their formated version
    col_mapping = dict(zip(col_names, col_names_formated))

    return col_mapping

def filter_data(data: pl.DataFrame, 
                col_keep: list, 
                year_stop: int, 
                flow_keep: list, 
                excluded_iso: list) -> pl.DataFrame:
    '''
    Function that filter uncomtrade data: select a set of columns for 
    readibility; filter so that dataframe does not exceed a certain year; keep
    flow type considered; exclude certain iso codes ~ countries; remove trade
    flows where reporter == partner; drop duplicates.

    Parameters
    ----------
    data : polars dataframe
        The trade data to filter.
    col_keep : list of strings
        The list of column names to keep.
    year_stop : integer
        The final year of the trade data considered.
    flow_keep : list of strings
        The trade flows to consider.
    excluded_iso : list of strings
        The list of ISO code to exclude from data.

    Returns
    -------
    filtered_data : polars dataframe
        The filtered data.
    '''

    filtered_data = (
        data
        .select(col_keep)
        # Format column names for homogenity
        .rename(format_col_names(col_names=col_keep))
        .filter(
            # Keep data in specified time range
            pl.col('period') <= str(year_stop-2),

            # Keep trade flow of interest
            pl.col('flow_code').is_in(flow_keep),

            # Remove non-country reporters and partners
            (~pl.col('reporter_iso')
            .str.contains('|'.join(excluded_iso))),
            (~pl.col('partner_iso')
            .str.contains('|'.join(excluded_iso))),

            # Delete "auto-" imports or exports (reporter = partner)
            pl.col('reporter_desc') != pl.col('partner_desc'),
        )
        # Drop potential duplicates
        .unique()
    )

    return filtered_data

def drop_outliers(data: pl.DataFrame,
                  percentiles: float = 0.05) :
    '''
    Function that drop outliers, i.e., values under the considered percentiles
    for net weight and primary value.

    Parameters
    ----------
    data : polars dataframe
        The trade data to filter.
    percentiles : float
        The percentiles to drop.

    Returns
    -------
    outliers_droped_data : polars dataframe
        The filtered data.
    '''
    stats_desc = (
        data
        .select(['net_wgt', 'primary_value'])
        .describe(percentiles=percentiles)
    )

    min_weight, min_value = (
        stats_desc
        .filter(pl.col('statistic').str.contains('%'))
        .select(['net_wgt', 'primary_value'])
    )

    outliers_droped_data = (
        data.filter(
            # Drop trade flow with net weight (kg) under fifth percentile
            ((pl.col('net_wgt') > min_weight.item()) | 
            (pl.col('net_wgt').is_null())),

            # Drop trade flow with value (USD) under fifth percentile
            ((pl.col('primary_value') > min_value.item()) | 
            (pl.col('primary_value').is_null()))
        )
    )

    return outliers_droped_data

# Log file edition
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

# Load data
due_diligence_data = pl.read_parquet(snakemake.input[0])
logging.info(f"Uncomtrade data for product under due diligence downloaded:\n {due_diligence_data}\n")

placebo_data = pl.read_parquet(snakemake.input[1])
logging.info(f"Placebo data downloaded:\n {placebo_data}\n")

due_diligence_codes_hs6 = pl.read_csv(snakemake.input[2], infer_schema=False)

placebo_codes_hs6 = pl.read_csv(snakemake.input[3], infer_schema=False)


logging.info("\n\n\n*** Processing due diligence data ***\n")

# Filter data
due_diligence_input = filter_data(
    data=due_diligence_data,
    col_keep=snakemake.params['col_keep'],
    year_stop=snakemake.params['year_stop'],
    flow_keep=snakemake.params['flow_keep'],
    excluded_iso=snakemake.params['excluded_iso']
)
logging.info(f"Filtered data:\n {due_diligence_input}\n")
logging.info(f"Filtered data stats:\n {due_diligence_input.describe()}\n")

# Drop outliers = values under fifth percentile for weight (kg) and value (USD)
due_diligence_input = drop_outliers(data=due_diligence_input)
logging.info(f"Saved data after deleting outliers:\n {due_diligence_input}")
logging.info(f"Saved data after deleting outliers stats:\n {due_diligence_input.describe()}")

# Join parent_code column
cols = due_diligence_input.columns
insert = [i+1 for i in range(len(cols)) if cols[i] == 'cmd_code'][0]
cols[insert:insert] = ['parent_code']

due_diligence_input = (
    due_diligence_input
    .join(
        due_diligence_codes_hs6.select(['code', 'parent_code']).unique(),
        left_on='cmd_code',
        right_on='code',
        how='left'
    )
    .select(cols)
)

# Add policy column
due_diligence_input = (
    due_diligence_input
    .join(
        (
            due_diligence_codes_hs6
            .select(['reporter_iso', 'year_implementation', 'code'])
            .with_columns(policy=pl.lit(1))
        ),
        left_on=['reporter_iso', 'period', 'cmd_code'],
        right_on=['reporter_iso', 'year_implementation', 'code'],
        how='left'
    )
    .with_columns(
        pl.col('policy').fill_null(0)
    )
)

# Sort data before saving
due_diligence_input = due_diligence_input.sort(['period', 
                                                'reporter_iso', 
                                                'partner_iso', 
                                                'cmd_code', 
                                                'primary_value'])
logging.info(f"Saved data after join:\n {due_diligence_input}")
logging.info(f"Saved data columns:\n {due_diligence_input.columns}")
logging.info(f"Saved data desc stats:\n {due_diligence_input.describe()}")

# Save input data
due_diligence_input.write_csv(
    snakemake.output[0],
    separator=','
    )


logging.info("\n\n\n*** Processing placebo data ***\n")

# Filter data
placebo_input = filter_data(
    data=placebo_data,
    col_keep=snakemake.params['col_keep'],
    year_stop=snakemake.params['year_stop'],
    flow_keep=snakemake.params['flow_keep'],
    excluded_iso=snakemake.params['excluded_iso']
)
logging.info(f"Filtered data:\n {placebo_input}\n")
logging.info(f"Filtered data stats:\n {placebo_input.describe()}\n")

# Drop outliers = values under fifth percentile for weight (kg) and value (USD)
placebo_input = drop_outliers(data=placebo_input)
logging.info(f"Saved data after deleting outliers:\n {placebo_input}")
logging.info(f"Saved data after deleting outliers stats:\n {placebo_input.describe()}")

# Join parent_code column
cols = placebo_input.columns
insert = [i+1 for i in range(len(cols)) if cols[i] == 'cmd_code'][0]
cols[insert:insert] = ['parent_code']

placebo_input = (
    placebo_input
    .join(
        placebo_codes_hs6.select(['code', 'parent_code']).unique(),
        left_on='cmd_code',
        right_on='code',
        how='left'
    )
    .select(cols)
)

# Add policy column
placebo_input = (
    placebo_input
    .join(
        (
            placebo_codes_hs6
            .select(['reporter_iso', 'year_implementation', 'code'])
            .with_columns(policy=pl.lit(1))
        ),
        left_on=['reporter_iso', 'period', 'cmd_code'],
        right_on=['reporter_iso', 'year_implementation', 'code'],
        how='left'
    )
    .with_columns(
        pl.col('policy').fill_null(0)
    )
)

# Sort data before saving
placebo_input = placebo_input.sort(['period', 
                                    'reporter_iso', 
                                    'partner_iso', 
                                    'cmd_code', 
                                    'primary_value'])
logging.info(f"Saved data after join:\n {placebo_input}")
logging.info(f"Saved data columns:\n {placebo_input.columns}")
logging.info(f"Saved data desc stats:\n {placebo_input.describe()}")

# Save input data
placebo_input.write_csv(
    snakemake.output[1],
    separator=','
    )
