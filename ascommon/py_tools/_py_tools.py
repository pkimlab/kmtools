import os
import tempfile
import logging


def unique(l):
    """Return a list of unique elements of `l`, preserving order.

    Examples
    --------
    >>> list(unique([1, 2, 3, 2, 1]))
    [1, 2, 3]
    """
    seen = set()
    yield from (x for x in l if x not in seen and not seen.add(x))


def configure_logging(
        level='info',
        format='[%(levelname)s] - [%(name)s]: %(message)s'):
    """Get a logger with basic configurations.

    .. note::

        This function is obsolete.
        It was created when I did not know how logging works.
    """
    level_dict = {
        'debug': logging.DEBUG,
        'info': logging.INFO,
        'warning': logging.WARNING,
        'error': logging.ERROR
    }
    # "set up the default handler" using basicConfig
    logging.basicConfig(level=level_dict[level], format=format)
    logging.getLogger().setLevel(level_dict[level])
    logging.getLogger('sqlalchemy.engine').setLevel(logging.WARNING)
    logging.getLogger("paramiko").setLevel(logging.WARNING)
    logging.debug('Done configuring logging!')


def set_temp_dir(tempdir=None):
    if tempdir is None:
        if 'TMPDIR' in os.environ:
            tempdir = os.environ['TMPDIR']
        else:
            raise ValueError("""\
You must either provide a 'tempdir' parameter or set the '$TMPDIR' environment variable!\
""")
    os.makedirs(tempdir, exist_ok=True)
    tempfile.tempdir = tempdir


def get_temp_dir():
    """Set the temporary directory to be used by `tempfile`.

    Examples
    --------
    >>> set_temp_dir('/tmp/_test')
    >>> get_temp_dir()
    '/tmp/_test'
    """
    return tempfile.gettempdir()
