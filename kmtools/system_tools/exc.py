import subprocess


# Database
class ArchiveError(Exception):
    """Error working with 7zip archive file."""

    def __init__(self, result, error_message, return_code):
        super(ArchiveError, self).__init__(result)
        self.error_message = error_message
        self.return_code = return_code


class ArchiveNotFoundError(ArchiveError):
    """7zip archive file not found."""

    pass


class SubprocessError(subprocess.SubprocessError):
    def __init__(self, command, host, stdout, stderr, returncode):
        self.command = command
        self.host = host
        self.stdout = stdout
        self.stderr = stderr
        self.returncode = returncode
