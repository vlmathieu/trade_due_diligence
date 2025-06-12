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
uncomtrade_data = pl.read_parquet('/Users/valentinmathieu/Desktop/wd/trade_due_diligence/resources/public/uncomtrade_data.parquet.gzip')
input_data = pl.read_parquet('/Users/valentinmathieu/Desktop/wd/trade_due_diligence/results/input/input_data.parquet.gzip')
input_data_eu = pl.read_parquet('/Users/valentinmathieu/Desktop/wd/trade_due_diligence/results/intermediary/input_data_eu.parquet.gzip')