import functools
import logging
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


def retry_database(fn):
    """Decorator to keep probing the database untill you succeed."""
    _check_exception = functools.partial(check_exception, valid_exc=sa.exc.OperationalError)
    r = retry(
        retry_on_exception=_check_exception, wait_exponential_multiplier=1000,
        wait_exponential_max=60000, stop_max_attempt_number=7)
    return r(fn)


def retry_subprocess(fn):
    _check_exception = functools.partial(check_exception, valid_exc=subprocess.SubprocessError)
    r = retry(
        retry_on_exception=_check_exception, wait_exponential_multiplier=1000,
        wait_exponential_max=60000, stop_max_attempt_number=7)
    return r(fn)


def retry_archive(fn):
    """Decorator to keep probing the database untill you succeed."""
    _check_exception = functools.partial(check_exception, valid_exc=system_tools.exc.ArchiveError)
    r = retry(retry_on_exception=_check_exception, wait_fixed=2000, stop_max_attempt_number=2)
    return r(fn)
