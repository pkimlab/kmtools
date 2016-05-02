import os
from importlib import reload

version_suffix = '_v9'
code_path = os.path.dirname(__file__)
data_path = os.path.join(code_path, 'data/')

db_cache_filename = data_path + 'db_cache.h5'


def configure_logging(
        level='info',
        format='%(levelname)s:%(name)s:%(message)s'):
    """Get a logger with basic configurations.
    """
    import logging
    reload(logging)
    level_dict = {
        'debug': logging.DEBUG,
        'info': logging.INFO,
        'warning': logging.WARNING,
        'error': logging.ERROR
    }
    logging.basicConfig(level=level_dict[level], format=format)
    logging.getLogger('sqlalchemy.engine').setLevel(logging.WARNING)
    logging.getLogger("paramiko").setLevel(logging.WARNING)
    logging.debug('Done configuring logging!')
