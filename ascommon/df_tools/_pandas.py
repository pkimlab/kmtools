import logging
from collections import Counter
import numpy as np


logger = logging.getLogger(__name__)


def split_df(df, n_chunks):
    """Split DataFrame `df` into `nchunks` chunks."""
    chunk_size = int(np.ceil(df.shape[0] / n_chunks))
    assert n_chunks * chunk_size >= df.shape[0]
    chunks = []
    for i in range(0, df.shape[0], chunk_size):
        chunks.append(df[i:i + chunk_size])
    assert len(chunks) == n_chunks
    return chunks


def remove_duplicate_columns(df, keep_first=True, add_suffix=False):
    """Remove duplicate columns.

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame([[1,2,3], [4,5,6], [7,8,9]], columns=['a', 'b', 'a'])

    >>> remove_duplicate_columns(df)
       a  b
    0  1  2
    1  4  5
    2  7  8
    >>> remove_duplicate_columns(df, keep_first=False)
       b  a
    0  2  3
    1  5  6
    2  8  9
    >>> remove_duplicate_columns(df, add_suffix=True)
       a  b a_2
    0  1  2   3
    1  4  5   6
    2  7  8   9
    """
    seen = Counter()
    keep_i = []
    keep_name = []

    def my_enumerate(columns):
        if keep_first:
            yield from enumerate(columns)
        else:
            yield from enumerate(reversed(columns))

    for i, column in my_enumerate(df.columns):
        if column not in seen:
            seen.update([column])
            keep_i.append(i)
            keep_name.append(column)
        else:
            if add_suffix:
                column_orig = column
                suffix_i = 1
                while column in seen:
                    suffix_i += 1
                    column = '{}_{}'.format(column_orig, suffix_i)
                logger.info("Renamed column '{}' to '{}'.".format(column_orig, column))
                keep_i.append(i)
                keep_name.append(column)
            else:
                logger.info("Removed column '{}' at position {}.".format(column, i))

    if not keep_first:
        keep_i = list(reversed([df.shape[1] - i - 1 for i in keep_i]))
        keep_name = list(reversed(keep_name))

    df = df.iloc[:, keep_i]
    df.columns = keep_name
    return df
