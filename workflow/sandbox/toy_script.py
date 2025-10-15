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
