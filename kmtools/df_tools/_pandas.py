import logging
import string
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


def get_reverse_column(column):
    """Change column 'xxx_1' to 'xxx_2' and vice versa.

    Examples
    --------
    >>> get_reverse_column('hello')
    'hello'
    >>> get_reverse_column('good_bye_1')
    'good_bye_2'
    >>> get_reverse_column('world_2')
    'world_1'
    """
    for suffix, other in [('_1', '_2'), ('_2', '_1')]:
        if column.endswith(suffix):
            return column[:-len(suffix)] + other
    return column


def reverse_columns(columns):
    """Flip the suffix of columns which end in '_1' and '_2'.

    Examples
    --------
    >>> reverse_columns(['a_1', 'a_2', 'b_2', 'b_1', 'c'])
    ['a_2', 'a_1', 'b_1', 'b_2', 'c']
    """
    return [get_reverse_column(c) for c in columns]


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
                logger.info("Renamed column '%s' to '%s'.", column_orig, column)
                keep_i.append(i)
                keep_name.append(column)
            else:
                logger.info("Removed column '%s' at position %s.", column, i)

    if not keep_first:
        keep_i = list(reversed([df.shape[1] - i - 1 for i in keep_i]))
        keep_name = list(reversed(keep_name))

    df = df.iloc[:, keep_i]
    df.columns = keep_name
    return df


def random_df(n_rows: int=100000) -> pd.DataFrame:
    """Return a random class:`pandas.DataFrame`.

    Can help with tests where you need a random DataFrame.
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
