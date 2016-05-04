import os
import hashlib
import fcntl
import re
from contextlib import contextmanager
import pandas as pd


DB_CACHE_FILENAME = None


def parse_connection_str(connection_str):
    match = (
        re.findall(
            '(\w*)://(\w*):(\w*)@(\d*\.\d*\.\d*\.\d*):?(\d*)/(\w*)',
            connection_str)
    )
    if len(match):
        db, username, password, ip, port, db_table = match[0]
    else:
        db, username, password, ip, port, db_table = [None, None, None, None, None, None]
    return db, username, password, ip, port, db_table


@contextmanager
def lock_tables(tablenames, engine):
    """Lock a list of tables so that other processes can't access them while you do your thing."""
    if type(tablenames) not in {list, tuple}:
        tablenames = [tablenames]
    try:
        engine.execute('set innodb_lock_wait_timeout=14400')
        engine.execute('lock tables ' + ' '.join([t + ' write' for t in tablenames]))
        yield
    except:
        raise
    finally:
        engine.execute('unlock tables')


@contextmanager
def open_exclusively(filename, mode='a'):
    try:
        fd = os.open(filename, os.O_CREAT | os.O_RDWR)
        fcntl.lockf(fd, fcntl.LOCK_EX)
        f = os.fdopen(fd, mode)
        yield f
    except:
        raise
    finally:
        f.close()


def _validate_db_cache_filename(db_cache_filename):
    if db_cache_filename is None:
        db_cache_filename = DB_CACHE_FILENAME
    if db_cache_filename is None:
        error_message = """\
You must provide `db_cache_filename` or set `DB_CACHE_FILENAME`.\
"""
        raise Exception(error_message)


@contextmanager
def open_hdf5_exclusively(filename, mode='a', db_cache_filename=None):
    _validate_db_cache_filename(db_cache_filename)
    try:
        fd = os.open(filename, os.O_CREAT | os.O_RDWR)
        fcntl.lockf(fd, fcntl.LOCK_EX)
        f = pd.HDFStore(DB_CACHE_FILENAME, mode)
        yield f
    except:
        raise
    finally:
        f.close()


class QueryNotInCache(Exception):
    pass


def _set_db_cache(table_name, table_data, db_cache_filename=None):
    _validate_db_cache_filename(db_cache_filename)
    with open_hdf5_exclusively(DB_CACHE_FILENAME, 'r+') as hdf:
        hdf[table_name] = table_data


def _get_db_cache(table_name, db_cache_filename=None):
    _validate_db_cache_filename(db_cache_filename)
    with open_hdf5_exclusively(DB_CACHE_FILENAME, 'r') as hdf:
        if table_name not in hdf:
            table_data = None
        else:
            table_data = hdf[table_name]
    return table_data


def read_sql_cached(sql_query, engine=None, db_cache_filename=None):
    _validate_db_cache_filename(db_cache_filename)
    hashed_sql_query = hashlib.sha224(sql_query).hexdigest()
    with open_hdf5_exclusively(DB_CACHE_FILENAME, 'r+') as hdf:
        if hashed_sql_query in hdf:
            result_df = hdf[hashed_sql_query]
        elif engine is not None:
            result_df = pd.read_sql(sql_query, engine)
            hdf[hashed_sql_query] = result_df
        else:
            result_df = None
    return result_df
