from snakemake.script import snakemake
import logging
import polars as pl

# Log file edition
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

# Build corresponding table

corr_table = pl.from_dicts(
    [

        {
            'cmdCode': 4403,
            'cmdDesc': 'Wood in the rough, whether or not stripped of bark or sapwood, or roughly squared.',
            'faoCode': 1861,
            'faoDesc': 'Roundwood'
        },

        {
            'cmdCode': 4407,
            'cmdDesc': 'Wood sawn or chipped lengthwise, sliced or peeled, whether or not planed, sanded or end-jointed, of a thickness exceeding 6 mm.',
            'faoCode': 1872,
            'faoDesc': 'Sawnwood'
        },

        {
            'cmdCode': 4412,
            'cmdDesc': 'Plywood, veneered panels and similar laminated wood.',
            'faoCode': 1640,
            'faoDesc': 'Plywood and LVL' # LVL = Laminated veneer lumber
        }
    ]
)

logging.info(f"Corresponding table built:\n {corr_table}")

# Save corresponding table
corr_table.write_csv(
    file= snakemake.output[0],
    separator=';'
)
logging.info("Data saved.")
