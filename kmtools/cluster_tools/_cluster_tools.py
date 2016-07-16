import re


def parse_connection_string(connection_string):
    """Split `connection_string` into database parameters.

    Examples
    --------
    >>> from pprint import pprint
    >>> pprint(parse_connection_string('mysql://root:@localhost'))
    {'db_name': '',
     'db_type': 'mysql',
     'host_ip': 'localhost',
     'host_port': '',
     'password': '',
     'socket': '',
     'username': 'root'}
    >>> pprint(parse_connection_string('mysql://root:root_pass@192.168.0.1:3306/test'))
    {'db_name': 'test',
     'db_type': 'mysql',
     'host_ip': '192.168.0.1',
     'host_port': '3306',
     'password': 'root_pass',
     'socket': '',
     'username': 'root'}
    """
    db_params = {}
    (db_params['db_type'], db_params['username'], db_params['password'],
     db_params['host_ip'], db_params['host_port'], db_params['db_name'],
     db_params['socket']) = (
        re.match(
            '^(\w*)://(\w*):(\w*)@(localhost|[0-9\.]*)(|:[0-9]*)(|\/\w*)(|\?unix_socket=.*)$',
            connection_string)
        .groups()
    )
    db_params['host_port'] = db_params['host_port'].strip(':')
    db_params['db_name'] = db_params['db_name'].strip('/')
    db_params['socket'] = db_params['socket'].partition('?unix_socket=')[-1]
    return db_params


def iterate_parameters(parameter_grid, _params=None):
    """.

    Parameters
    ----------
    parameter_grid : dict
        Keys are parameters. Values are lists of values that the parameters can take.

    Examples
    --------
    >>> from pprint import pprint
    >>> pprint(list(iterate_parameters({'a': [1, 2], 'b': [3, 4]})))
    [{'a': 1, 'b': 3}, {'a': 1, 'b': 4}, {'a': 2, 'b': 3}, {'a': 2, 'b': 4}]
    """
    param_grid = parameter_grid.copy()
    # Don't modify dictionaries in-place
    if _params is None:
        _params = dict()
    # Terminal case
    if not param_grid:
        yield _params
        return
    # Recurse
    key, values = param_grid.popitem()
    for value in values:
        _params[key] = value
        try:
            yield from iterate_parameters(param_grid.copy(), _params.copy())
        except StopIteration:
            continue
