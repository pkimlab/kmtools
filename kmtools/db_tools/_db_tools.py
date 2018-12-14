import fcntl
import hashlib
import os
import re
from contextlib import contextmanager
from textwrap import dedent
from typing import NamedTuple, Optional, Union

import pandas as pd

DB_CACHE_FILENAME = None


class ConOpts(NamedTuple):
    """
    .. deprecated::
        Use `urllib.parse.ParseResult` instead.
    """

    name: str
    username: str
    password: Optional[str]
    url: str
    port: Optional[int]
    schema: str
    socket: str


def parse_connection_string(connection_string):
    """Split `connection_string` into a dictionary of connection properties.

    Notes:
        The returned dictionary maps everything to strings (including `port`)
        to make it compatible with :class:`configparser.configParser`.

    Args:
        connection_string: String describing database connection in SQLAlchemy-compatible format.

    .. deprecated::
        Use `urllib.parse.urlparse` instead.

    Returns:
        A tuple of connection parameters.
    """
    db_params = {}
    (
        db_params["name"],
        db_params["username"],
        db_params["password"],
        db_params["url"],
        db_params["port"],
        db_params["schema"],
        db_params["socket"],
    ) = re.match(
        r"^(\w*)"  # name
        r"://"
        r"(|\w*)"  # username
        r"(|:\w*)"  # password
        r"(|@localhost|@[a-zA-Z0-9\.-]*|@[0-9\.]*)"  # url
        r"(|:[0-9]*)"  # port
        r"(|\/[^?]*)"  # schema
        r"(|\?unix_socket=.*)$",  # socket
        connection_string,
    ).groups()
    if db_params["password"].startswith(":"):
        if db_params["password"] == ":":
            db_params["password"] = ""
        else:
            db_params["password"] = db_params["password"][1:]
    else:
        db_params["password"] = None
    db_params["url"] = db_params["url"].lstrip("@")
    db_params["port"] = db_params["port"].lstrip(":")
    try:
        db_params["port"] = int(db_params["port"])
    except ValueError:
        db_params["port"] = None
    db_params["schema"] = (
        db_params["schema"][1:] if db_params["schema"].startswith("/") else db_params["schema"]
    )
    db_params["socket"] = db_params["socket"].partition("?unix_socket=")[-1]
    return ConOpts._make((db_params[key] for key in ConOpts._fields))


def make_connection_string(vargs: Union[dict, ConOpts]):
    """Join a dictionary of connection properties (`vargs`) into a connection string.

    .. deprecated::
        Use `urllib.parse.urlunparse` / `urllib.parse.urlunsplit` instead.
    """
    if isinstance(vargs, ConOpts):
        vargs = vargs._asdict()
    vargs["password"] = (
        ":{}".format(vargs["password"])
        if vargs.get("password") is not None and not vargs.get("schema", "").startswith("/")
        else ""
    )
    vargs["url"] = "@{}".format(vargs["url"]) if vargs.get("url") else ""
    vargs["port"] = ":{}".format(vargs["port"]) if vargs.get("port") else ""
    vargs["schema"] = "/{}".format(vargs["schema"]) if vargs.get("schema") else "/"
    vargs["socket"] = "?unix_socket={}".format(vargs["socket"]) if vargs.get("socket") else ""
    connection_string = "{name}://{username}{password}{url}{port}{schema}{socket}".format(**vargs)
    return connection_string


@contextmanager
def lock_tables(tablenames, engine):
    """Lock a list of tables so that other processes can't access them while you do your thing."""
    if type(tablenames) not in {list, tuple}:
        tablenames = [tablenames]
    try:
        engine.execute("set innodb_lock_wait_timeout=14400")
        engine.execute("lock tables " + " ".join([t + " write" for t in tablenames]))
        yield
    except Exception:
        raise
    finally:
        engine.execute("unlock tables")


def _validate_db_cache_filename(db_cache_filename):
    if db_cache_filename is None:
        db_cache_filename = DB_CACHE_FILENAME
    if db_cache_filename is None:
        error_message = dedent(
            """\
            You must provide `db_cache_filename` or set `DB_CACHE_FILENAME`.
            """
        ).strip()
        raise Exception(error_message)


@contextmanager
def open_hdf5_exclusively(filename, mode="a", db_cache_filename=None):
    _validate_db_cache_filename(db_cache_filename)
    try:
        fd = os.open(filename, os.O_CREAT | os.O_RDWR)
        fcntl.lockf(fd, fcntl.LOCK_EX)
        f = pd.HDFStore(DB_CACHE_FILENAME, mode)
        yield f
    except Exception:
        raise
    finally:
        f.close()


class QueryNotInCache(Exception):
    pass


def _set_db_cache(table_name, table_data, db_cache_filename=None):
    _validate_db_cache_filename(db_cache_filename)
    with open_hdf5_exclusively(DB_CACHE_FILENAME, "r+") as hdf:
        hdf[table_name] = table_data


def _get_db_cache(table_name, db_cache_filename=None):
    _validate_db_cache_filename(db_cache_filename)
    with open_hdf5_exclusively(DB_CACHE_FILENAME, "r") as hdf:
        if table_name not in hdf:
            table_data = None
        else:
            table_data = hdf[table_name]
    return table_data


def read_sql_cached(sql_query, engine=None, db_cache_filename=None):
    _validate_db_cache_filename(db_cache_filename)
    hashed_sql_query = hashlib.sha224(sql_query).hexdigest()
    with open_hdf5_exclusively(DB_CACHE_FILENAME, "r+") as hdf:
        if hashed_sql_query in hdf:
            result_df = hdf[hashed_sql_query]
        elif engine is not None:
            result_df = pd.read_sql(sql_query, engine)
            hdf[hashed_sql_query] = result_df
        else:
            result_df = None
    return result_df
