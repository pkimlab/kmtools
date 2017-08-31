import contextlib
import logging

import paramiko

logger = logging.getLogger(__name__)


@contextlib.contextmanager
def ssh_client(ip):
    try:
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(ip)
        yield ssh
    except Exception as e:
        logger.error(type(e))
        logger.error(e)


class SSHClient:

    def __init__(self, ssh_host):
        self.ssh_host = ssh_host
        self.ssh = paramiko.SSHClient()
        self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    def __enter__(self):
        logger.debug("Initializing SSH client: '{}'".format(self.ssh_host))
        hostname, username, password = self._parse_ssh_connection()
        self.ssh.connect(hostname, username=username, password=password)
        return self.ssh

    def _parse_ssh_connection(self):
        hostname = self.ssh_host
        username = None
        password = None
        if '@' in hostname:
            username, hostname = self.ssh_host.split('@')
            if ':' in username:
                username, password = username.split(':')
        return hostname, username, password

    def __exit__(self, exc_type, exc_value, exc_tb):
        self.ssh.close()
        if exc_type or exc_value or exc_tb:
            import traceback
            logger.error(exc_type)
            logger.error(exc_value)
            logger.error(traceback.print_tb(exc_tb))
            return True
        else:
            return False
