"""Opinionated way of loading or dumping dataframe(s) as DSV files.

.. note::

    You should probably use `feather.write_dataframe(df)` instead.
"""
import csv
import logging
import os
import os.path as op

import numpy as np
import pandas as pd

logger = logging.getLogger()


def dump(what, where, how=".tsv.gz"):
    """Dump one or more :class:`pandas.DataFrame` objects to one or more files.

    Parameters
    ----------
    what: pd.DataFrame | dict[pd.DataFrame]
        What should be saved?
    where: str
        Where should it be saved (i.e. path)?
    how: str, default='.tsv.gz'
        How should it be saved (i.e. what extension to use)?
    """
    if isinstance(what, pd.DataFrame):
        return _dump_single(what, where, how)
    elif isinstance(what, dict) and all(isinstance(v, pd.DataFrame) for v in what.values()):
        return _dump_multiple(what, where, how)
    else:
        raise NotImplementedError("dump({}, {}, {how})")


def _dump_single(what, where, how):
    file = where + how
    sep = _guess_sep(file)
    compression = _guess_compression(file)
    what.to_csv(
        file,
        sep=sep,
        index=False,
        compression=compression,
        na_rep="\\N",
        quoting=csv.QUOTE_NONNUMERIC,
    )
    what.dtypes.to_pickle(file + ".dtype")


def _dump_multiple(what, where, how):
    for key, df in what.items():
        file = op.join(where, key)
        dump(df, file, how)


def load(where, how=".tsv.gz"):
    """Load one or more :class:`pandas.DataFrame` objects from one or more files.

    Parameters
    ----------
    where: str
        Where should the files be loaded from?
    how: str, default='.tsv.gz'
        What is the format of the file (i.e. what extension does it have)?
    """
    if op.isfile(where + how):
        return _load_single(where, how)
    elif op.isdir(where):
        return _load_multiple(where, how)
    else:
        raise Exception("File or folder '{}' does not exist!".format(where))


def _load_single(where, how):
    file = where + how
    sep = _guess_sep(file)
    compression = _guess_compression(file)
    dtypes = pd.read_pickle(file + ".dtype")
    dtypes.loc[dtypes == "<M8[ns]"] = np.dtype("O")
    df = pd.read_csv(
        file,
        sep=sep,
        compression=compression,
        dtype=dtypes.to_dict(),
        na_values=["\\N"],
        na_filter=False,
        quoting=csv.QUOTE_NONNUMERIC,
    )
    return df


def _load_multiple(where, how):
    dfs = {}
    for filename in os.listdir(where):
        if not filename.endswith(how):
            continue
        filename = filename.rpartition(how)[0]
        dfs[filename] = _load_single(op.join(where, filename), how)
    return dfs


def _guess_sep(file):
    if ".tsv" in file:
        return "\t"
    else:
        return ","


def _guess_compression(filename):
    if filename.endswith(".gz"):
        return "gzip"
    elif filename.endswith(".bz2"):
        return "bz2"
    elif filename.endswith(".xz"):
        return "xz"
    else:
        return None
