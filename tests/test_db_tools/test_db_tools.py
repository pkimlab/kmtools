import pytest

from kmtools.db_tools import ConOpts, make_connection_string, parse_connection_string

TEST_DATA = [
    # (connection_string, ConOpts)
    ('sqlite:////absolute/path/to/foo.db', ConOpts('sqlite', '', None, '', None,
                                                   '/absolute/path/to/foo.db', '')),
    ('mysql://user:@localhost/', ConOpts('mysql', 'user', '', 'localhost', None, '', '')),
    ('mysql://user:pass@192.168.0.1:3306/test', ConOpts('mysql', 'user', 'pass', '192.168.0.1',
                                                        3306, 'test', '')),
    ('mysql://user@192.168.0.1:3306/test?unix_socket=/tmp/mysql.sock', ConOpts(
        'mysql', 'user', None, '192.168.0.1', 3306, 'test', '/tmp/mysql.sock')),
]


@pytest.mark.parametrize("connection_string_, con_opts_", TEST_DATA)
def test_connection_string_ops(connection_string_, con_opts_):
    con_opts = parse_connection_string(connection_string_)
    assert con_opts == con_opts_
    connection_string = make_connection_string(con_opts)
    assert connection_string == connection_string_
