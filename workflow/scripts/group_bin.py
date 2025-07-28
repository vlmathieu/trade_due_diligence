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

# Get due diligence policies and correspondong countries from parameters
pol_params = snakemake.params['policies']
policies = list(pol_params.keys())
logging.info(f"Due diligence policies considered: {policies}.\n")
policy_countries = [_ for p in policies for _ in pol_params[p]['countries']]
logging.info(f"Countries that implemented due diligence policies: {policy_countries}\n")

# Get list of risky countries
risky = snakemake.params['risky']
logging.info(f"List of risky countries: {risky}\n")

# Produce binary variables for policies
for _ in policies:
    input_data = (
        input_data
        .with_columns(pl.col('period').cast(pl.Int32))
        .with_columns(
                pl.when(pl.col('period') >= pol_params[_]['year'])
                .then(1)
                .otherwise(0)
                .alias(f'{_}')
            )
    )

# Produce binary variables for country groups
for _ in ['reporter', 'partner']:
    input_data = (
        input_data
        .with_columns(
            pl.when(pl.col(f'{_}_desc').is_in(risky))
            .then(1)
            .otherwise(0)
            .alias(f'{_}_risky')
        )
        .with_columns(
            pl.when(pl.col(f'{_}_desc').is_in(policy_countries))
            .then(1)
            .otherwise(0)
            .alias(f'{_}_policy')
        )
        .with_columns(
            pl.when(~pl.col(f'{_}_desc').is_in(policy_countries+risky))
            .then(1)
            .otherwise(0)
            .alias(f'{_}_non_risky')
        )
    )