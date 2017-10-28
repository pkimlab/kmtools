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
    parquet_file = op.join(basedir, '{}.parquet'.format(table_name))
    table_pq = pa.Table.from_pandas(table_df.rename(columns=format_column))
    parquet_schema = table_pq.schema
    parquet_writer = pq.ParquetWriter(parquet_file, parquet_schema, compression='snappy')
    parquet_writer.write_table(table_pq)
    return parquet_file


def convert_hdf5_to_parquet(h5_file, parquet_file, chunksize=100000):

    stream = pd.read_csv(h5_file, chunksize=chunksize)

    for i, chunk in enumerate(stream):
        print("Chunk {}".format(i))

        if i == 0:
            # Infer schema and open parquet file on first chunk
            parquet_schema = pa.Table.from_pandas(df=chunk).schema
            parquet_writer = pq.ParquetWriter(parquet_file, parquet_schema, compression='snappy')

        table = pa.Table.from_pandas(chunk, schema=parquet_schema)
        parquet_writer.write_table(table)

    parquet_writer.close()
