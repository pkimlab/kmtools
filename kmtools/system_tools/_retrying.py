import functools
import logging
import socket
import subprocess

import sqlalchemy as sa
from retrying import retry

from kmtools import system_tools

logger = logging.getLogger(__name__)


def check_exception(exc, valid_exc):
    logger.error('The following exception occured:\n{}'.format(exc))
    to_retry = isinstance(exc, valid_exc)
    if to_retry:
        logger.error('Retrying...')
    return to_retry


def _retry_urlopen(fn):
    """Retry downloading data from a url after a timeout.

    Examples
    --------
    >>> import urllib.request
    >>> with urllib.request.urlopen('http://google.ca') as ifh:
    ...     data = _retry_urlopen(ifh.read)()
    """
    _check_exception = functools.partial(check_exception, valid_exc=socket.timeout)
    wrapper = retry(
        retry_on_exception=_check_exception,
        wait_exponential_multiplier=1000,
        wait_exponential_max=10000,
        stop_max_attempt_number=5)
    return wrapper(fn)


def retry_database(fn):
    """Decorator to keep probing the database untill you succeed."""
    _check_exception = functools.partial(check_exception, valid_exc=sa.exc.OperationalError)
    wrapper = retry(
        retry_on_exception=_check_exception,
        wait_exponential_multiplier=1000,
        wait_exponential_max=60000,
        stop_max_attempt_number=7)
    return wrapper(fn)


def retry_subprocess(fn):
    _check_exception = functools.partial(check_exception, valid_exc=subprocess.SubprocessError)
    wrapper = retry(
        retry_on_exception=_check_exception,
        wait_exponential_multiplier=1000,
        wait_exponential_max=60000,
        stop_max_attempt_number=7)
    return wrapper(fn)


def retry_archive(fn):
    """Decorator to keep probing the database untill you succeed."""
    _check_exception = functools.partial(check_exception, valid_exc=system_tools.exc.ArchiveError)
    wrapper = retry(
        retry_on_exception=_check_exception, wait_fixed=2000, stop_max_attempt_number=2)
    return wrapper(fn)
