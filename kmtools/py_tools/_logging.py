"""Enable '{}' formatting for logging statements.

Source: `Use of alternative formatting styles`_.

_`Use of alternative formatting styles`:
    https://docs.python.org/3/howto/logging-cookbook.html#use-of-alternative-formatting-styles
"""
import functools
import logging

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
