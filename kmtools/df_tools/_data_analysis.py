"""Data Analysis Tools.
"""
import fcntl
import logging
import os
from contextlib import contextmanager

import numpy as np
import pandas as pd
import tables

logger = logging.getLogger(__name__)


def print2(a, b, *args, x=60):
    template = "{:%d}{:,}" % x
    formatted_template = template.format(a, b)
    for arg in args:
        formatted_template += " " + str(arg)
    print(formatted_template)


def print_full(x):
    """Print the entire Dataframe / Series."""
    pd.set_option("display.max_columns", None)
    pd.set_option("display.max_rows", None)
    print(x)
    pd.reset_option("display.max_columns")
    pd.reset_option("display.max_rows")


def remove_duplicates(seq, keep_null=True):
    if keep_null:
        seen = set()
    else:
        seen = set([np.nan, None, ""])
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


@contextmanager
def open_hdf5_exclusively(filename, mode="a"):
    try:
        fd = os.open(filename, os.O_CREAT | os.O_RDWR)
        fcntl.lockf(fd, fcntl.LOCK_EX)
        f = tables.openFile(filename, mode)
        yield f
    except Exception:
        raise
    finally:
        f.close()
