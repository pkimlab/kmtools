import os.path as op

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from . import format_column


def write_parquet_tables(basedir='', **kwargs):
    files = []
    for table_name, table_df in kwargs.popitem():
        assert isinstance(table_df, pd.DataFrame)
        files.append(_write_parquet_table(table_name, table_df, basedir=basedir))
    return files


def _write_parquet_table(table_name, table_df, basedir=''):
    parquet_file = op.join(basedir, f'{table_name}.parquet')
    table_pq = pa.Table.from_pandas(table_df.rename(columns=format_column))
    parquet_schema = table_pq.schema
    parquet_writer = pq.ParquetWriter(parquet_file, parquet_schema, compression='snappy')
    parquet_writer.write_table(table_pq)
    return parquet_file
