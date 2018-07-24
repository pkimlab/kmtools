import functools
import logging
import os
import sys
import threading
import time
from contextlib import contextmanager

getLogger = logging.getLogger

LOGGING_LEVELS = {
    0: logging.ERROR,
    1: logging.WARNING,
    2: logging.INFO,
    3: logging.DEBUG,
}


class LogPipe(threading.Thread):
    """Redirect messages written to a pipe into a log.

    Source: https://codereview.stackexchange.com/a/17959/62763

    Examples:
        >>> import subprocess
        >>> from contextlib import closing
        >>> with closing(LogPipe(logging.DEBUG)) as log_pipe:
        ...     cp = subprocess.run(
        ...         ["echo", "hello world"], stdout=log_pipe, universal_newlines=True)
    """

    def __init__(self, fn):
        """Setup the object with a logger and a loglevel and start the thread.
        """
        super().__init__()
        self.daemon = False
        self._fn = fn
        self._fout, self._fin = os.pipe()
        self._pipe_reader = os.fdopen(self._fout)
        self.start()

    def fileno(self):
        """Return the write file descriptor of the pipe.
        """
        return self._fin

    def run(self):
        """Run the thread, logging everything.
        """
        for line in iter(self._pipe_reader.readline, ''):
            line = line.strip()
            if line:
                self._fn(line)
        self._pipe_reader.close()

    def close(self):
        """Close the write end of the pipe.
        """
        # Close the input channel
        os.close(self._fin)
        # Wait for the output channel to flush out
        while self.is_alive():
            time.sleep(0.01)


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
    """Enable '{}' formatting for logging statements.

    Source: `Use of alternative formatting styles`_.

    _`Use of alternative formatting styles`:
        https://docs.python.org/3/howto/logging-cookbook.html#use-of-alternative-formatting-styles
    """

    def __init__(self, logger, extra=None):
        super(StyleAdapter, self).__init__(logger, extra or {})

    def log(self, level, msg, *args, **kwargs):
        if self.isEnabledFor(level):
            msg, kwargs = self.process(msg, kwargs)
            self.logger._log(level, Message(msg, args), (), **kwargs)


class LoggingContext(object):

    def __init__(self, logger, level=None, handler=None, close=True):
        self.logger = logger
        self.level = level
        self.handler = handler
        self.close = close

    def __enter__(self):
        if self.level is not None:
            self.old_level = self.logger.level
            self.logger.setLevel(self.level)
        if self.handler:
            self.logger.addHandler(self.handler)

    def __exit__(self, et, ev, tb):
        if self.level is not None:
            self.logger.setLevel(self.old_level)
        if self.handler:
            self.logger.removeHandler(self.handler)
        if self.handler and self.close:
            self.handler.close()
        # implicit return of None => don't swallow exceptions


class WritableObject:
    """A writable object which writes everything to the logger."""

    def __init__(self, logger):
        self.logger = logger

    def write(self, string):
        self.logger.debug(string.strip())


@functools.wraps(logging.getLogger)
def get_logger(*args, **kwargs):
    return StyleAdapter(getLogger(*args, **kwargs))


def patch_getLogger():
    logging.getLogger = get_logger


def log_function_calls(logger):
    """Log every call of the decorated function."""

    def decorator(fn):

        @functools.wraps(fn)
        def wrapper(*args, **kwargs):
            logger.warning(fn.__name__ + '(' + ', '.join(args) +
                           ', '.join('{}={}'.format(k, v) for k, v in kwargs.items()) + ')')
            return fn(*args, **kwargs)

        return wrapper

    return decorator


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
    except Exception:
        raise
    finally:
        sys.stdout = original_stdout
        for i in range(len(logger.handlers)):
            logger.handlers[i].formatter = original_formatters[i]
