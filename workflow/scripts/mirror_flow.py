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
                flow_keep: list = ['M', 'X'], 
                excluded_iso: list = ['XX', '_X']) -> pl.DataFrame:
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
        The trade flows to consider. Default value: ['M', 'X']
    excluded_iso : list of strings
        The list of ISO code to exclude from data. Default value: ['XX', '_X']

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

def one_side_flows(data: pl.dataframe.frame.DataFrame, 
                   flow_code: str,
                   flow_type: list = ['wgt', 'value'],
                   trader_role: list = ['reporter', 'partner']):
    '''
    Function that extracts and formats one part of the mirror flows (either
    exporter reports or importer reports) and return one-side trade flows.

    Parameters
    ----------
    data : polars dataframe
        The UNComtrade data.
    flow_code : string
        Code that specifies if a trade flow report is an import or an export. 
        Must be 'M' for import or 'X' for export.
    flow_type : list of string
        Type of flows considered for building mirror flows. The items of 
        flow_type must be contained in column names. Default value is 
        ['wgt', 'value'] for commodities weight and trade values.
    trader_role : list of string
        The role of reporter or partner of traders in trade flow reporting.
        Default value is ['reporter', 'partner'].
    
    Raises
    ------
    ValueError
        Raises ValueError is flowDesc is not 'M' or 'X'.

    Returns
    -------
    extraction_formated : polars dataframe
        Either import- or export-side trade flows formated to build mirro flows.

    '''
    
    # Values of flow_code must be 'M' for import or 'X' for export, error is not
    valid = {'M','X'}
    if flow_code not in valid:
        raise ValueError('one_side_flow: flow_code must be one of %r' % valid)
    
    # List column names of flow types
    flows = [_ for _ in data.columns if any(s in _ for s in flow_type)]
    traders = [_ for _ in data.columns if any(s in _ for s in trader_role)]

    # Format one_side_flow columns names depending on flow_code
    if flow_code == 'X':
        mapping = (
            # If export, reporter is the exporter and partner the importer
            dict(zip(traders,
                     [_.replace('reporter', 'exporter')
                       .replace('partner', 'importer') for _ in traders])) |
            # Add suffix _exp to flow type columns
            dict(zip(flows, 
                     [f'{_}_exp' for _ in flows]))
        )
    elif flow_code == 'M':
        mapping = (
            # If import, reporter is the importer and partner the exporter
            dict(zip(traders,
                     [_.replace('reporter', 'importer')
                       .replace('partner', 'exporter') for _ in traders])) |
            # Add suffix _imp to flow type columns
            dict(zip(flows, 
                     [f'{_}_imp' for _ in flows]))
        )
    # Extract data and format columns
    extraction_formated = (
        data
        .filter(pl.col('flow_code') == flow_code)
        .rename(mapping)
        .drop('flow_code')
    )
    
    return extraction_formated

def compute_mirror_flows(data: pl.dataframe.frame.DataFrame,
                         report_suffix: list = ['_imp', '_exp']):
    '''
    Function that assembles one-side trade flows into mirror flows based on 
    cleaned UN Comtrade data.

    Parameters
    ----------
    data : polars dataframe
        The UNComtrade data.
    report_suffix : list of string
        List of suffixes contained in column names referring to country reports. 
        Default is ['_imp', '_exp'] for importer or exporter reports, 
        respectively.

    Returns
    -------
    mirror_flows : polars dataframe
        The trade mirror flows.

    '''
    
    exp_flow = one_side_flows(data,'X')
    imp_flow = one_side_flows(data,'M')
    join_col = [_ for _ in exp_flow.columns 
                if not any(s in _ for s in report_suffix)]

    mirror_flows = exp_flow.join(imp_flow,
                                 on=join_col,
                                 how='full',
                                 coalesce=True)
    
    return mirror_flows

# Log file edition
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
    
# Load input trade data
due_diligence_data = pl.read_parquet(snakemake.input[0])

# Create mirror flows
due_diligence_mirror_flows = compute_mirror_flows(
    filter_data(
        due_diligence_data,
        col_keep=snakemake.params['col_keep'], 
        year_stop=snakemake.params['year_stop'], 
        excluded_iso= snakemake.params['excluded_iso']
    )
)
logging.info(f"\nMirror flows:\n {due_diligence_mirror_flows}\n")
logging.info(f"\nMirror flows description:\n {due_diligence_mirror_flows.describe()}\n")

# Save mirror flows
due_diligence_mirror_flows.write_csv(
    snakemake.output[0],
    separator=','
    )
