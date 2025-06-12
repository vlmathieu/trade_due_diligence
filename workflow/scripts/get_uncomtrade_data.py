from snakemake.script import snakemake
import polars as pl
import polars.selectors as cs
import comtradeapicall
from tqdm import tqdm 

def get_uncomtrade(apikey: str, year: str, cmd: str, flow: str):
    '''
    Function that downloads UN Comtrade data for a given year, a given 
    commodity, and a given trade flow. Need an API key.

    Parameters
    ----------
    apikey : sting
        The API subscription key to download data.
    year : string
        The year of trade.
    cmd : string
        The commodity code.
    flow : string
        The trade flow to download (import, export, re-import, re-export...).

    Returns
    -------
    data : polars dataframe
        The UN Comtrade data for a given year, commodity, and trade flow.

    '''
    
    uncomtrade_data = comtradeapicall.getFinalData(
        apikey,
        typeCode        = 'C',          # typeCode(str) : Product type. Goods (C) or Services (S)
        freqCode        = 'A',          # freqCode(str) : The time interval at which observations occur. Annual (A) or Monthly (M)
        clCode          = 'HS',         # clCode(str) : Indicates the product classification used and which version (HS, SITC)
        period          = year,         # period(str) : Combination of year and month (for monthly), year for (annual)
        reporterCode    = None,         # reporterCode(str) : The country or geographic area to which the measured statistical phenomenon relates
        cmdCode         = cmd,          # cmdCode(str) : Product code in conjunction with classification code
        flowCode        = flow,         # flowCode(str) : Trade flow or sub-flow (exports, re-exports, imports, re-imports, etc.)
        partnerCode     = None,         # partnerCode(str) : The primary partner country or geographic area for the respective trade flow
        partner2Code    = None,         # partner2Code(str) : A secondary partner country or geographic area for the respective trade flow
        customsCode     = None,         # customsCode(str) : Customs or statistical procedure
        motCode         = None,         # motCode(str) : The mode of transport used when goods enter or leave the economic territory of a country
        format_output   = 'JSON',       # format_output(str) : The output format. CSV or JSON
        breakdownMode   = 'classic',    # breakdownMode(str) : Option to select the classic (trade by partner/product) or plus (extended breakdown) mode
        includeDesc     = True          # includeDesc(bool) : Option to include the description or not
        )
        
    return uncomtrade_data

def chunks(lst: list, n: int):
    '''
    Yield successive n-sized chunks from a list.

    Parameters
    ----------
    lst : list
        List to divide into chunks.
    n : integer
        Size of the chunks.
    '''
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def get_uncomtrade_bulk(apikey: str, years: list, cmdCode: list, flowCode: list):
    '''
    Function that downloads UN Comtrade data for a several years, commodities, 
    and trade flows. Need an API key.

    Parameters
    ----------
    apikey : string
        The API subscription key to download data.
    years : list of strings
        The years of trade.
    cmdCode : list of strings
        The commodity codes.
    flowCode : list of strings
        The trade flow to download (import, export, re-import, re-export...).

    Returns
    -------
    uncomtrade_data : polars dataframe
        The UN Comtrade data for a several years, commodity, and trade flows.

    '''

    # Download data by batch of 5 years and 3 commodities
    uncomtrade_years_batch = (
        [pl.from_pandas(
            get_uncomtrade(
                str(apikey),
                ','.join([str(_) for _ in years_batch]),
                ','.join([str(_) for _ in cmd_batch]),
                ','.join([str(_) for _ in flowCode])
            ))
        for years_batch in chunks(years, 5)
        for cmd_batch in tqdm(chunks(cmdCode, 3))]
    )

    # Concatenate all batch
    uncomtrade_data = pl.concat(
        [df for df in uncomtrade_years_batch if df.shape != (0,0)],
        how='vertical_relaxed'
    )

    # Check of years, commodities, and different flows considered
    check_list = (
        sorted(set(uncomtrade_data['period'].unique())) == sorted(set([str(_) for _ in years])),
        sorted(set(uncomtrade_data['cmdCode'].unique())) == sorted(set([str(_) for _ in cmdCode])),
        sorted(set(uncomtrade_data['flowCode'].unique())) == sorted(set([str(_) for _ in flowCode]))
    )

    return uncomtrade_data, check_list

UN_Comtrade_data, check_list = get_uncomtrade_bulk(
    snakemake.params['apikey'],
    list(range(snakemake.params['year_start'], snakemake.params['year_stop'])),
    snakemake.params['cmdCode'],
    snakemake.params['flowCode']
)

print("\nDataframe head: \n\n", UN_Comtrade_data.head(5), "\n")
print("\nDataframe size (rows, columns): ", UN_Comtrade_data.shape, "\n")

# Save data if check list passed
if all(check_list):
    print('Data have been checked.\n')   
    UN_Comtrade_data.write_parquet(
        snakemake.output[0],
        compression='gzip'
        )
else:
    print('Issues found in data download.\n')
