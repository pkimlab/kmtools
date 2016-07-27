import os
import hashlib
import fcntl
import re
from contextlib import contextmanager
import pandas as pd


DB_CACHE_FILENAME = None


def parse_connection_string(connection_string):
    """Split `connection_string` into a dictionary of connection properties.

    Examples
    --------
    >>> from pprint import pprint
    >>> pprint(parse_connection_string('mysql://user:@localhost'))
    {'db_password': '',
     'db_port': '',
     'db_schema': '',
     'db_socket': '',
     'db_type': 'mysql',
     'db_url': 'localhost',
     'db_username': 'user'}
    >>> pprint(parse_connection_string('mysql://user:pass@192.168.0.1:3306/test'))
    {'db_password': 'pass',
     'db_port': '3306',
     'db_schema': 'test',
     'db_socket': '',
     'db_type': 'mysql',
     'db_url': '192.168.0.1',
     'db_username': 'user'}
    >>> pprint(parse_connection_string('sqlite:////absolute/path/to/foo.db'))
    {'db_password': '',
     'db_port': '',
     'db_schema': '/absolute/path/to/foo.db',
     'db_socket': '',
     'db_type': 'sqlite',
     'db_url': '',
     'db_username': ''}
    >>> connection_string = 'mysql://user:pass@192.168.0.1:3306/test?unix_socket=/tmp/mysql.sock'
    >>> pprint(parse_connection_string(connection_string))
    {'db_password': 'pass',
     'db_port': '3306',
     'db_schema': 'test',
     'db_socket': '/tmp/mysql.sock',
     'db_type': 'mysql',
     'db_url': '192.168.0.1',
     'db_username': 'user'}
    """
    db_params = {}
    (db_params['db_type'], db_params['db_username'], db_params['db_password'],
     db_params['db_url'], db_params['db_port'], db_params['db_schema'],
     db_params['db_socket']) = (
        re.match(
            '^(\w*)'  # db_type
            '://'
            '(|\w*:)'  # db_username
            '(|\w*)'  # db_password
            '(|@localhost|@[0-9\.]*)'  # db_url
            '(|:[0-9]*)'  # db_port
            '(|\/[^?]*)'  # db_schema
            '(|\?unix_socket=.*)$',  # db_socket
            connection_string)
        .groups()
    )
    db_params['db_username'] = db_params['db_username'].rstrip(':')
    db_params['db_url'] = db_params['db_url'].lstrip('@')
    db_params['db_port'] = db_params['db_port'].lstrip(':')
    db_params['db_schema'] = (
        db_params['db_schema'][1:]
        if db_params['db_schema'].startswith('/')
        else db_params['db_schema'])
    db_params['db_socket'] = db_params['db_socket'].partition('?unix_socket=')[-1]
    return db_params


def make_connection_string(**vargs):
    """Join a dictionary of connection properties (`vargs`) into a connection string.

    Examples
    --------
    >>> make_connection_string(**{ \
        'db_password': '', \
        'db_port': '', \
        'db_schema': '', \
        'db_socket': '', \
        'db_type': 'mysql', \
        'db_url': 'localhost', \
        'db_username': 'user'})
    'mysql://user:@localhost/'
    >>> make_connection_string(**{ \
        'db_password': 'pass', \
        'db_port': '3306', \
        'db_schema': 'test', \
        'db_socket': '', \
        'db_type': 'mysql', \
        'db_url': '192.168.0.1', \
        'db_username': 'user'})
    'mysql://user:pass@192.168.0.1:3306/test'
    >>> make_connection_string(**{ \
        'db_password': '', \
        'db_port': '', \
        'db_schema': '/absolute/path/to/foo.db', \
        'db_socket': '', \
        'db_type': 'sqlite', \
        'db_url': '', \
        'db_username': ''})
    'sqlite:////absolute/path/to/foo.db'
    """
    if vargs['db_username']:
        vargs['db_username'] = vargs['db_username'] + ':'
    if vargs['db_url']:
        vargs['db_url'] = '@' + vargs['db_url']
    if vargs['db_port']:
        assert vargs['db_url']
        vargs['db_port'] = ':' + vargs['db_port']
    if vargs['db_schema'] is not None:
        vargs['db_schema'] = '/' + vargs['db_schema']
    if vargs['db_socket']:
        vargs['db_socket'] = '?unix_socket=' + vargs['db_socket']
    connection_string = (
        '{db_type}://{db_username}{db_password}{db_url}{db_port}{db_schema}{db_socket}'
        .format(**vargs)
    )
    return connection_string


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
