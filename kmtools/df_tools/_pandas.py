import string
import logging
from collections import Counter, OrderedDict
import numpy as np
import pandas as pd


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


# Save and load CSV files with one line
def dump_csv(df: pd.DataFrame, file: str) -> None:
    """Export DataFrame as a CSV file with dtype annotations.
    """
    sep = _guess_sep(file)
    compression = _guess_compression(file)
    df.to_csv(file, sep=sep, index=False, compression=compression)
    df.dtypes.to_pickle(file + '.dtype')


def load_csv(file: str) -> pd.DataFrame:
    """Import DataFrame from a CSV file and dtype annotations.
    """
    sep = _guess_sep(file)
    compression = _guess_compression(file)
    dtypes = pd.read_pickle(file + '.dtype')
    df = pd.read_csv(file, sep=sep, compression=compression, dtype=dtypes.to_dict())
    return df


def _guess_sep(file):
    if '.tsv' in file:
        return '\t'
    else:
        return ','


def _guess_compression(file):
    if file.endswith('.gz'):
        return 'gzip'
    elif file.endswith('.bz2'):
        return 'bz2'
    elif file.endswith('.xz'):
        return 'xz'
    else:
        return None


# Helps with tests
def random_df(n_rows: int=100000, n_cols: int=5) -> pd.DataFrame:
    """Return a random class:`pandas.DataFrame`.
    """
    fns = [
        lambda: np.array(
            [''.join(np.random.choice(list(string.printable)) for _ in range(12))] * n_rows
        ),
        lambda: np.array(
            [''.join(np.random.choice(list(string.printable)) for _ in range(256))] * n_rows
        ),
        lambda: np.array(
            [''.join(np.random.choice(list(string.printable)) for _ in range(1024))] * n_rows
        ),
        lambda: np.random.randint(0, 10000000000, n_rows),
        lambda: np.random.randn(n_rows),
    ]
    cols = OrderedDict()
    for name, fn in zip(string.ascii_uppercase, fns):
        cols[name] = fn()
    return pd.DataFrame(cols)
