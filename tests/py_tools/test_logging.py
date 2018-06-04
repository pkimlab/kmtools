import logging
import os
import subprocess
from contextlib import closing
from unittest.mock import patch

import pytest

from kmtools import py_tools


@pytest.mark.parametrize("level", [logging.DEBUG, logging.INFO, logging.WARNING])
def test_log_pipe(level):
    py_tools.LogPipe.logger.setLevel(logging.INFO)

    msg = 'test simple write'
    with patch('kmtools.py_tools.LogPipe.logger') as mock_logger:
        with closing(py_tools.LogPipe(level)) as log_pipe:
            os.write(log_pipe.fileno(), msg.encode())
    if level >= logging.INFO:
        mock_logger.log.assert_called_once_with(level, msg)
    else:
        mock_logger.lock.assert_not_called()

    msg = 'test subprocess'
    with patch('kmtools.py_tools.LogPipe.logger') as mock_logger:
        with closing(py_tools.LogPipe(level)) as log_pipe:
            subprocess.run(['echo', msg], stdout=log_pipe, universal_newlines=True)
    if level >= logging.INFO:
        mock_logger.log.assert_called_once_with(level, msg)
    else:
        mock_logger.lock.assert_not_called()
