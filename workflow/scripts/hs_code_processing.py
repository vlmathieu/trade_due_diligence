from snakemake.script import snakemake
import logging
import polars as pl

def hs_disaggregate(codes, hs, level_up: str):
    '''
    Function that disaggregates HS codes from aggregated levels (HS2 or HS4) to
    disaggregated level (HS4 or HS6).

    Parameters
    ----------
    codes : polars dataframe
        The dataframe containing the HS codes to disaggregate.
    hs : polars dataframe
        The dataframe of all HS codes at all levels for all classifications.
    level_up : string
        The level to disaggregate.

    Returns
    -------
    codes_low : polars dataframe
        A dataframe with the disaggregated HS codes.
    '''

    # Columns to keep
    col_keep = ['country',
                'due_diligence_policy',
                'year_signature',
                'year_implementation',
                'classification',
                'code']

    # Get all HS codes at the upper level
    parent_codes = (codes.filter(pl.col('level') == level_up)
                         .select(pl.col('code'))
                         .unique()
                         .to_series()
                         .to_list())

    # Get corresponding disaggregated HS codes (one level lower)
    codes_up_to_low = (hs.filter(pl.col('parent_code').is_in(parent_codes))
                         .select(pl.exclude('description')))

    # Filter HS codes to keep row at the upper level
    codes_up = (codes.filter(pl.col('level') == level_up)
                     .select(col_keep)
                     .rename({'code': 'parent_code'}))

    # Join with the disaggregated HS codes
    codes_low = codes_up.join(
        codes_up_to_low,
        left_on=['classification', 'parent_code'],
        right_on=['classification', 'parent_code']
    )
    
    # Reorder columns
    codes_low = codes_low.select(codes.columns)

    return codes_low

def country_details(codes, desc: list, iso: list):
    '''
    Function that disaggregates HS codes country attribute, from countries 
    unions- such as the European Union- to individual countries and join the 
    countries' name and iso.

    Parameters
    ----------
    codes : polars dataframe
        The dataframe containing the HS codes to disaggregate.
    desc : list of string
        The list of all country names, for countries of interest.
    level_up : string
        The list of all country ISO codes, for countries of interest.

    Returns
    -------
    codes_low : polars dataframe
        A dataframe with the disaggregated HS codes.
    '''

    # Create dataframe of countries details
    countries = pl.DataFrame(
        {
            'reporter_desc': desc,
            'reporter_iso': iso
        }
    )

    # Dataframe of due diligence policy per country group
    policy_group = (
        codes
        .select(['country', 'due_diligence_policy'])
        .unique()
    )

    # Dataframe of due diligence policy per individual country
    policy_ind = (
        countries
        .join(
            policy_group,
            left_on='reporter_desc',
            right_on='country',
            how='left'
        )
        .fill_null('European Union Timber Regulation')
    )

    # Disaggregated HS codes for individual countries
    codes_low = (
        codes
        .join(
            policy_ind,
            on='due_diligence_policy'
        )
        .select(pl.exclude('country'))
    )
    
    # Reorder columns
    cols = ['reporter_iso', 'reporter_desc'] + codes.columns[1:-1]
    codes_low = codes_low.select(cols)

    return codes_low

# Log file edition
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

# Load data
hs_classification = pl.read_csv(
    snakemake.input[0],
    infer_schema=False
)
due_diligence_codes = pl.read_csv(
    snakemake.input[1],
    infer_schema=False
)
placebo_codes = pl.read_csv(
    snakemake.input[2],
    infer_schema=False
).filter(pl.col('parent_code').is_in(snakemake.params['placebo_codes']))
logging.info(f"\nHS classification:\n {hs_classification}\n")
logging.info(f"\nDue diligence HS codes:\n {due_diligence_codes}\n")
logging.info(f"\nPlacebo HS codes:\n {placebo_codes}\n")


logging.info("\n\n\n*** Processing due diligence data ***\n")

# Disaggregate due diligence HS4 codes to HS6 codes
due_diligence_hs4_to_hs6 = hs_disaggregate(
    codes=due_diligence_codes,
    hs=hs_classification,
    level_up='4'
)
logging.info(f"\nDue diligence HS4 to HS6 codes: \n {due_diligence_hs4_to_hs6}\n")

# Concatenate with due diligence HS6 codes
due_diligence_hs6 = pl.concat(
    [
        due_diligence_codes.filter(pl.col('level') == '6'),
        due_diligence_hs4_to_hs6,
    ],
    how="vertical"
).sort(['year_signature', 'year_implementation', 'code'])
logging.info(f"\nDue diligence HS6 codes: \n {due_diligence_hs6}\n")

logging.info(f"\nInitial dataframe size (rows, columns): {due_diligence_codes.shape}\n")
logging.info(f"\nNew dataframe size (rows, columns): {due_diligence_hs6.shape}\n")

# Disaggregate due diligence HS6 codes at the individual country level
due_diligence_hs6 = country_details(
    codes=due_diligence_hs6,
    desc=snakemake.params['desc'],
    iso=snakemake.params['iso']
)
logging.info(f"\nDue diligence HS6 codes at individual country level: \n {due_diligence_hs6}\n")
logging.info(f"\nList individual countries: \n {sorted(due_diligence_hs6.select('reporter_desc').unique().to_series().to_list())}\n")

# Save processed data
due_diligence_hs6.write_csv(
    snakemake.output[0],
    separator=','
)


logging.info("\n\n\n*** Processing placebo data ***\n")

# Disaggregate placebo HS4 codes to HS6 codes
placebo_hs6 = hs_disaggregate(
    codes=placebo_codes,
    hs=hs_classification,
    level_up='4'
).sort(['year_signature', 'year_implementation', 'code'])
logging.info(f"\nPlacebo HS6 codes: \n {placebo_hs6}\n")

# Disaggregate placebo HS6 codes at the individual country level
placebo_hs6 = country_details(
    codes=placebo_hs6,
    desc=snakemake.params['desc'],
    iso=snakemake.params['iso']
)
logging.info(f"\nPlacebo HS6 codes at individual country level: \n {placebo_hs6}\n")
logging.info(f"\nList individual countries: \n {sorted(placebo_hs6.select('reporter_desc').unique().to_series().to_list())}\n")

# Save processed data
placebo_hs6.write_csv(
    snakemake.output[1],
    separator=','
)
