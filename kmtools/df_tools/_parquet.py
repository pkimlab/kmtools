import inspect
import logging
from contextlib import closing
from pathlib import Path
from typing import Any, Callable, Dict, Generator, Union

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

logger = logging.getLogger(__name__)


def csv_to_parquet(
    csv_file: str,
    parquet_file: str,
    csv_kwargs: Dict[str, Any] = None,
    parquet_kwargs: Dict[str, Any] = None,
):
    """Convert a CSV file to an Apache Parquet file."""
    dfs = pd.read_csv(csv_file, **(csv_kwargs or {}))
    write_parquet_file(dfs, parquet_file, **(parquet_kwargs or {}))


def write_parquet_file(
    dfs: Union[pd.DataFrame, Generator[pd.DataFrame, None, None]],
    parquet_file: Union[str, Path],
    **kwargs,
) -> None:
    """Write a DataFrame, or an iterator over DataFrames, to a Parquet file.

    Args:
        dfs: A DataFrame, or an iterator over DataFrames, that are to be written.
        parquet_file: The name of the Parquet file to which to write the DataFrame(s).
        **kwargs: Keyword arguments that will get passed to `pa.Table.from_pandas`,
            `pq.ParquetWriter`, and `pq.ParquetWriter.write_table` methods.

    Raises:
        Execption: If keyword arguments don't match the signatures of the methods
            listed above.
    """
    if isinstance(parquet_file, Path):
        parquet_file = parquet_file.absolute().as_posix()

    kwargs.setdefault("flavor", "spark")
    kwargs.setdefault("compression", "snappy")

    pa_kwargs = {k: kwargs.pop(k) for k in list(kwargs) if k in ["preserve_index"]}
    # This should work in the future versions of PyArrow:
    # pw_kwargs = extract_kwargs(pq.ParquetWriter, kwargs)
    pw_kwargs = {
        k: kwargs.pop(k)
        for k in list(kwargs)
        if k
        in ["flavor", "version", "use_dictionary", "compression", "use_deprecated_int96_timestamps"]
    }
    pq_kwargs = extract_kwargs(pq.ParquetWriter.write_table, kwargs)
    if kwargs:
        raise Exception(
            "Not all arguments where used during the call to _get_parser! "
            f"Remaining kwargs: {kwargs}"
        )

    if isinstance(dfs, pd.DataFrame):
        table = pa.Table.from_pandas(dfs, **pa_kwargs)
    else:
        table = pa.Table.from_pandas(next(dfs), **pa_kwargs)

    with closing(pq.ParquetWriter(parquet_file, table.schema, **pw_kwargs)) as parquet_writer:
        logger.debug("Writing data to file %s...", parquet_file)
        parquet_writer.write_table(table, **pq_kwargs)
        if not isinstance(dfs, pd.DataFrame):
            for i, chunk in enumerate(dfs):
                logger.debug("Writing chunk number %i...", i + 1)
                table = pa.Table.from_pandas(chunk, **pa_kwargs)
                parquet_writer.write_table(table, **pq_kwargs)


def extract_kwargs(fn: Callable[..., Any], kwargs: Dict[str, Any]) -> Dict[str, Any]:
    """Extract keyword arguments accepted by `fn` from `kwargs`.

    Notes:
        * This function mutates the `kwargs` argument
          by removing extracted keyword arguments.
    """
    fn_params = set(inspect.signature(fn).parameters)
    fn_kwargs = {k: kwargs.pop(k) for k in list(kwargs) if k in fn_params}
    return fn_kwargs
