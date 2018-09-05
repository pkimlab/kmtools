import functools
import logging
import os.path as op
import random
import string
import tempfile

# from hypothesis import given
import hypothesis.strategies as st
import numpy as np
import pandas as pd
import pytest
from pandas.util.testing import assert_frame_equal

import kmtools.df_tools

logger = logging.getLogger()

_how = [sep + compression for sep in (".csv", ".tsv") for compression in ("", ".gz", ".bz2", ".xz")]


@pytest.fixture(scope="module", params=_how)
def how(request):
    return request.param


def test_single(how):
    """Test different separators and compression methods for a single DSV file."""
    logger.debug("test_single(%s)", how)
    df = _gen_df()
    filename = _gen_filename()
    with tempfile.TemporaryDirectory() as tmpdir:
        file = op.join(tmpdir, filename)
        kmtools.df_tools.dump(df, file, how)
        df_out = kmtools.df_tools.load(file, how)
    logger.debug("df1:\n%s\ndf2:\n%s", df.head(), df_out.head())
    assert_frame_equal(df, df_out)


def test_multiple(how):
    """Test different separators and compression methods for multiple DSV files."""
    logger.debug("test_multiple(%s)", how)
    dfs = {_gen_filename(): _gen_df() for _ in range(12)}
    with tempfile.TemporaryDirectory() as tmpdir:
        kmtools.df_tools.dump(dfs, tmpdir, how)
        dfs_out = kmtools.df_tools.load(tmpdir, how)
    assert not set(dfs) ^ set(dfs_out)
    for key in dfs:
        logger.debug("df1:\n%s\ndf2:\n%s", dfs[key].head(), dfs_out[key].head())
        assert_frame_equal(dfs[key], dfs_out[key])


def _assert_allclose_object(df1, df2):
    logger.debug("_assert_allclose_object(%s, %s)", df1.head(), df2.head())
    assert (
        (
            df1.select_dtypes(include=[object]).fillna(0).values
            == df2.select_dtypes(include=[object]).fillna(0).values
        )
        .all()
        .all()
    )


def _assert_allclose_numeric(df1, df2):
    logger.debug("_assert_allclose_numberic(%s, %s)", df1.head(), df2.head())
    assert df1.select_dtypes(exclude=[object]).equals(df2.select_dtypes(exclude=[object]))


def _gen_df(n_rows=1000):
    """Generate a random DataFrame.

    .. note::

        Generating a random value for every row takes too long,
        so we use the same value for every row.
        This may not be the safest way to do it though!
    """
    logger.debug("_gen_df(%s)", n_rows)

    _get_finite_value = functools.partial(_get_numeric_value, filt=np.isfinite)
    data = {
        "integer": np.repeat(
            _get_finite_value(
                st.integers(
                    min_value=np.iinfo(np.int64).min * 0.99, max_value=np.iinfo(np.int64).max * 0.99
                )
            ),
            n_rows,
        ),
        "float": np.repeat(
            _get_finite_value(
                st.floats(
                    min_value=np.finfo(np.float64).min * 0.99,
                    max_value=np.finfo(np.float64).max * 0.99,
                )
            ),
            n_rows,
        ),
        "bool": np.repeat(_get_finite_value(st.booleans()), n_rows),
        "text": np.repeat(
            st.text(alphabet=string.ascii_letters, min_size=0, max_size=100).example(), n_rows
        ),
    }
    df = pd.DataFrame(data)
    # print(df.head())
    return df


def _get_numeric_value(hst, filt):
    _counter = 0
    while True:
        try:
            return hst.filter(filt).example()
        except TypeError:
            _counter += 1
            if _counter <= 100:
                continue
            raise


def _gen_filename(n_char=64):
    """Generate a random filename."""
    logger.debug("_gen_filename(%s)", n_char)
    return "".join([random.choice(string.ascii_letters) for _ in range(n_char)])
