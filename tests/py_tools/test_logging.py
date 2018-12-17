import os
import subprocess
from contextlib import closing
from unittest import mock

from kmtools import py_tools


def test_log_pipe_1():
    mock_logger = mock.MagicMock()
    mock_logger("starting")
    with closing(py_tools.LogPipe(mock_logger)) as log_pipe:
        os.write(log_pipe.fileno(), "test".encode())
    mock_logger("done")
    mock_logger.assert_has_calls([mock.call("starting"), mock.call("test"), mock.call("done")])


def test_log_pipe_2():
    mock_logger = mock.MagicMock()
    mock_logger("starting")
    with closing(py_tools.LogPipe(mock_logger)) as log_pipe:
        subprocess.run(["echo", "test"], stdout=log_pipe, universal_newlines=True)
    mock_logger("done")
    mock_logger.assert_has_calls([mock.call("starting"), mock.call("test"), mock.call("done")])
