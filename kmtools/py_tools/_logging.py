"""Enable '{}' formatting for logging statements.

Source: `Use of alternative formatting styles`_.

_`Use of alternative formatting styles`:
    https://docs.python.org/3/howto/logging-cookbook.html#use-of-alternative-formatting-styles
"""
import functools
import logging
import sys
from contextlib import contextmanager

getLogger = logging.getLogger


class Message(object):
    def __init__(self, fmt, args):
        self.fmt = fmt
        self.args = args

    def __str__(self):
        if isinstance(self.fmt, str):
            return self.fmt.format(*self.args)
        else:
            return str(self.fmt)


class StyleAdapter(logging.LoggerAdapter):
    def __init__(self, logger, extra=None):
        super(StyleAdapter, self).__init__(logger, extra or {})

    def log(self, level, msg, *args, **kwargs):
        if self.isEnabledFor(level):
            msg, kwargs = self.process(msg, kwargs)
            self.logger._log(level, Message(msg, args), (), **kwargs)


@functools.wraps(logging.getLogger)
def get_logger(*args, **kwargs):
    return StyleAdapter(getLogger(*args, **kwargs))


def patch_getLogger():
    logging.getLogger = get_logger


# ========= Log print statements =========

class WritableObject:
    """A writable object which writes everything to the logger."""

    def __init__(self, logger):
        self.logger = logger

    def write(self, string):
        self.logger.debug(string.strip())


@contextmanager
def log_print_statements(logger):
    """Channel print statements to the debug logger.

    Useful for modules that default to printing things instead of using a logger (Modeller...).
    """
    original_stdout = sys.stdout
    original_formatters = []
    for i in range(len(logger.handlers)):
        original_formatters.append(logger.handlers[0].formatter)
        logger.handlers[i].formatter = logging.Formatter('%(message)s')
    wo = WritableObject(logger)
    try:
        sys.stdout = wo
        yield
    except:
        raise
    finally:
        sys.stdout = original_stdout
        for i in range(len(logger.handlers)):
            logger.handlers[i].formatter = original_formatters[i]
