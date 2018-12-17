import functools
import logging
import os
import re
import shlex
import signal
import subprocess
import time
import warnings
from contextlib import contextmanager
from typing import Optional

import psutil

from ._remote import SSHClient
from ._retrying import retry_subprocess

logger = logging.getLogger(__name__)


def format_system_command(system_command):
    """Remove repeating whitespace from system command."""
    return " ".join(re.split(" *", system_command.strip()))


def _set_process_group(parent_process_group_id):
    """Set group_id of the child process to the group_id of the parent process.

    This way when you delete the parent process you also delete all the children.
    """
    child_process_id = os.getpid()
    os.setpgid(child_process_id, parent_process_group_id)


@functools.wraps(subprocess.run)
def run(system_command, cwd=None, **vargs):
    """Wrapper over `subprocess.run` with more reasonable defaults.

    Used by ELAPSAM among others, so do not remove...
    """
    logger.debug(system_command)
    if not isinstance(system_command, (list, tuple)):
        system_command = shlex.split(system_command)
    process = subprocess.run(
        system_command,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        preexec_fn=lambda: _set_process_group(os.getpgrp()),
        **vargs,
    )
    process.stdout = process.stdout.strip()
    process.stderr = process.stderr.strip()
    if process.returncode != 0:
        logger.info(process.stdout)
        logger.error(process.stderr)
        raise subprocess.CalledProcessError(
            process.returncode, system_command, process.stdout, process.stderr
        )
    return process


def start_subprocess(system_command):
    process = subprocess.Popen(
        shlex.split(system_command),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        bufsize=1,
    )
    return process


def iter_stdout(p):
    for line in p.stdout:
        line = line.strip()
        if " [Note] " in line:
            line = line.partition(" [Note] ")[-1]
        if not line:
            if p.poll() is None:
                time.sleep(0.1)
                continue
            else:
                # logger.debug("DONE! (reached an empty line)")
                return
        yield line


def execute(system_command: str, cwd: Optional[str] = None) -> None:
    """Execute a system command, passing STDERR to logger.

    Source: https://stackoverflow.com/a/4417735/2063031

    .. deprecated:: Use `kmtools.py_tools.LogPipe` instead.
    """
    warnings.warn(
        "`kmtools.system_tools.execute` is deprecated; use `kmtools.py_tools.LogPipe` instead.",
        DeprecationWarning,
    )
    logger.info("system_command: '%s'", system_command)
    popen = subprocess.Popen(
        shlex.split(system_command),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        cwd=cwd,
        universal_newlines=True,
    )
    for stdout_line in iter(popen.stdout.readline, ""):
        logger.debug(stdout_line.strip())
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, system_command)


@retry_subprocess
def run_command(system_command, host=None, allowed_returncodes=[0]):
    """Run system command either locally or over ssh.

    Notes:
        * For consistency with remote execution, local commands are executed with ``shell=True``.
          Accordingly, we should follow all neccessary precautions.
    """
    # === Run command ===
    logger.debug(system_command)
    if host is not None:
        stdout, stderr, returncode = _run_command_ssh(system_command, host)
    else:
        stdout, stderr, returncode = _run_command_local(system_command)
    # === Process results ===
    if returncode not in allowed_returncodes:
        raise subprocess.CalledProcessError(
            returncode=returncode, cmd=system_command, output=stdout, stderr=stderr
        )
    elif "warning" in stdout.lower() or "error" in stdout.lower():
        logger.warning("Command ran with warnings / errors! (%s)", stdout)
    else:
        logger.debug("Command ran successfully! (%s)", stdout)
    return stdout, stderr, returncode


def _run_command_ssh(system_command, host):
    logger.debug("Running on host: '%s'", host)
    with SSHClient(host) as ssh:
        _stdin, _stdout, _stderr = ssh.exec_command(system_command)
        stdout = _stdout.read().decode().strip()
        stderr = _stderr.read().decode().strip()
        returncode = _stdout.channel.recv_exit_status()
    return stdout, stderr, returncode


def _run_command_local(system_command):
    logger.debug("Running locally")
    process = subprocess.run(
        system_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True,
        universal_newlines=True,
    )
    return process.stdout.strip(), process.stderr.strip(), process.returncode


def kill_child_processes(parent_pid, sig=signal.SIGTERM):
    try:
        parent = psutil.Process(parent_pid)
    except psutil.NoSuchProcess:
        logger.debug("Counld not find parent process with pid '%s'", parent_pid)
        return
    children = parent.children(recursive=True)
    for process in children:
        logger.debug("Killing process '%s'...", process)
        process.send_signal(sig)


@contextmanager
def print_heartbeats():
    """Spawn a fork that prints a message every minute.

    This is required for travis-ci.
    """
    from elaspic import conf

    # Don't print random stuff if not testing
    if not conf.CONFIGS["testing"]:
        yield
        return
    # Print a heartbeat to keep travis happy.
    pid = os.fork()
    if pid == 0:
        while True:
            time.sleep(60)
            logger.info("Subprocess is still running...")
        os._exit()
    try:
        yield
    finally:
        os.kill(pid, 15)
        os.waitpid(pid, 0)
