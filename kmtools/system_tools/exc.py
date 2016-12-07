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
