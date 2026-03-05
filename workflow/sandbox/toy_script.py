import polars as pl
import polars.selectors as cs
import pickle
import re
import numpy as np
import networkx as nx
from scipy.stats import kurtosis
from scipy.stats import skew
import pandas as pd
import itertools

# Parameters
year_start = 2020
year_stop = 2024
cmdCode = [4403, 4407, 4408, 4812]
flowCode = ['M', 'X']
excluded_iso = ['XX', '_X', '\\d']

# Load data
uncomtrade_data = pl.read_parquet('/Users/valentinmathieu/Desktop/wd/trade_due_diligence/resources/public/uncomtrade.parquet.gzip')
input_data = pl.read_parquet('/Users/valentinmathieu/Desktop/wd/trade_due_diligence/results/input/input_uncomtrade.parquet.gzip')
input_data_eu = pl.read_parquet('/Users/valentinmathieu/Desktop/wd/trade_due_diligence/results/intermediary/input_uncomtrade_eu.parquet.gzip')

input_data.write_csv(
    '/Users/valentinmathieu/Desktop/input_uncomtrade.csv',
    separator=';'
)

dic = {'lacey_act': {'year': 2008, 'country': ['USA'], 'iso': ['USA']}, 'ilpa': {'year': 2012, 'country': ['Australia'], 'iso': ['AUS']}, 'eutr': {'year': 2013, 'country': ['Austria', 'Belgium', 'Belgium-Luxembourg (...1998)', 'Bulgaria', 'Croatia', 'Cyprus', 'Czechia', 'Denmark', 'Estonia', 'Finland', 'France', 'Germany', 'Greece', 'Hungary', 'Ireland', 'Italy', 'Latvia', 'Lithuania', 'Luxembourg', 'Malta', 'Netherlands', 'Poland', 'Portugal', 'Romania', 'Slovakia', 'Slovenia', 'Spain', 'Sweden', 'United Kingdom'], 'iso': ['AUT', 'BEL', 'BGR', 'CYP', 'CZE', 'DEU', 'DNK', 'ESP', 'EST', 'FIN', 'FRA', 'GBR', 'GRC', 'HRV', 'HUN', 'IRL', 'ITA', 'LTU', 'LUX', 'LVA', 'MLT', 'NLD', 'POL', 'PRT', 'ROU', 'SVK', 'SVN', 'SWE']}, 'cwa': {'year': 2016, 'country': ['Japan'], 'iso': ['JPN']}, 'asut': {'year': 2018, 'country': ['Rep. of Korea'], 'iso': ['KOR']}}
policies = list(dic.keys())
dic['lacey_act']['country']

for _ in policies:
    print(_)

p = 'lacey_act'
dic[p]
test = input_data
for _ in policies:
    test = (
        test
            .with_columns(
                pl.col('period').cast(pl.Int32)
            )
            .with_columns(
                pl.when(
                    (pl.col('period') >= dic[_]['year']) & 
                    (pl.col('reporter_desc').is_in(dic[_]['country'])))
                .then(1)
                .otherwise(0)
                .alias(f'reporter_{_}')
            )
            .with_columns(
                pl.when(
                    (pl.col('period') >= dic[_]['year']) & 
                    (pl.col('partner_desc').is_in(dic[p]['country'])))
                .then(1)
                .otherwise(0)
                .alias(f'partner_{_}')
            )
    )

###########

# Function
def hs_disaggregate(hs_codes, hs, level_up: str):
    '''
    Function that disaggregates HS codes from aggregated levels (HS2 or HS4) to
    disaggregated level (HS4 or HS6).

    Parameters
    ----------
    hs_codes : polars dataframe
        The dataframe containing the HS codes to disaggregate.
    hs : polars dataframe
        The dataframe of all HS codes at all levels for all classifications.
    level_up : string
        The level to disaggregate.

    Returns
    -------
    hs_codes_low : polars dataframe
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
    codes_up = (hs_codes.filter(pl.col('level') == level_up)
                .select(pl.col('code'))
                .unique()
                .to_series()
                .to_list())

    # Get corresponding disaggregated HS codes (one level lower)
    hs_codes_up_to_low = (hs.filter(pl.col('parent_code').is_in(codes_up))
                            .select(pl.exclude('description')))

    # Filter HS codes to keep row at the upper level
    hs_codes_up = (hs_codes.filter(pl.col('level') == level_up)
                   .select(col_keep))

    # Join with the disaggregated HS codes
    hs_codes_low = hs_codes_up.join(
        hs_codes_up_to_low,
        left_on=['classification', 'code'],
        right_on=['classification', 'parent_code']
    )
    
    return hs_codes_low

# Load data
hs = pl.read_csv(
    "resources/public/hs_code.csv",
    infer_schema=False
)
due_diligence_hs_raw = pl.read_csv(
    "resources/inhouse/due_diligence_hs_code_raw.csv",
    infer_schema=False
)

# Data processing

hs2_codes = (due_diligence_hs_raw
    .filter(pl.col('level') == '2')
    .select(pl.col('code'))
    .unique()
    .to_series()
    .to_list()
)

hs2_to_hs4_codes = (hs
    .filter(pl.col('parent_code').is_in(hs2_codes))
    .select(pl.col('code'))
    .unique()
    .to_series()
    .to_list()
)

hs4_codes = list(
    set(
        hs2_to_hs4_codes + 
        due_diligence_hs_raw
        .filter(pl.col('level') == '4')
        .select(pl.col('code'))
        .unique()
        .to_series()
        .to_list()
    )
)

hs4_to_hs6_codes = (hs
    .filter(pl.col('parent_code').is_in(hs4_codes))
    .select(pl.exclude('description'))
)

col_keep = ['country', 'due_diligence_policy', 'year_signature', 'year_implementation', 'classification',
            'code']
due_diligence_hs4 = (due_diligence_hs_raw
    .filter(pl.col('level') == '4')
    .select(col_keep)
)

test = due_diligence_hs4.join(
    hs4_to_hs6_codes,
    left_on=['classification', 'code'],
    right_on=['classification', 'parent_code']
)


hs_codes = due_diligence_hs_raw
hs = hs
level_up = '2'

# Columns to keep
col_keep = ['country',
            'due_diligence_policy',
            'year_signature',
            'year_implementation',
            'classification',
            'code']

# Get all HS codes at the upper level
hs_codes_up = (hs_codes.filter(pl.col('level') == level_up)
                .select(pl.col('code'))
                .unique()
                .to_series()
                .to_list())

# Get corresponding disaggregated HS codes (one level lower)
hs_codes_up_to_low = (hs.filter(pl.col('parent_code').is_in(hs_codes_up))
                        .select(pl.exclude('description')))

# Filter HS codes to keep row at the upper level
hs_codes_filter_up = (due_diligence_hs_raw.filter(pl.col('level') == level_up)
                      .select(col_keep))

# Join with the disaggregated HS codes
hs_codes_low = hs_codes_filter_up.join(
    hs_codes_up_to_low,
    left_on=['classification', 'code'],
    right_on=['classification', 'parent_code']
)

hs_codes_low

test = hs_disaggregate(hs_codes=due_diligence_hs_raw, hs=hs, level_up='2')