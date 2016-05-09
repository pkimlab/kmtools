from importlib import reload


CACHE_DIR = None  # place to store temporary files


def get_cache_dir(cache_dir=None):
    if cache_dir is None:
        cache_dir = CACHE_DIR
    if cache_dir is None:
        raise EnvironmentError('`ascommon.settings.CACHE_DIR` is not defined!')
    return cache_dir


def configure_logging(
        level='info',
        format='%(levelname)s:%(name)s:%(message)s'):
    """Get a logger with basic configurations."""
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
