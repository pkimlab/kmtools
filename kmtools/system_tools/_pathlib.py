import logging
import os
import os.path as op
import shutil
import string
from contextlib import contextmanager

logger = logging.getLogger(__name__)


# =============================================================================
# Filename string formatting
# =============================================================================

def slugify(filename_string):
    valid_chars = "-_.()" + string.ascii_letters + string.digits
    return ''.join(c if c in valid_chars else '_' for c in filename_string)


def remove_extensions(filename, extensions):
    """Remove extensions from file.

    Examples
    --------
    >>> remove_extensions('/tmp/a/b/c.d.e.f.g', ['.e', '.f', '.g'])
    '/tmp/a/b/c.d'
    >>> remove_extensions('do.re.mi', ['do', 'mi'])
    'do.re'
    """
    # Add missing '.'
    extensions = [(ext if ext.startswith('.') else '.' + ext) for ext in extensions]
    # Strip extensions
    while True:
        file_name, file_ext = op.splitext(filename)
        if file_ext in extensions:
            filename = file_name
        else:
            break
    return filename


def strip_ps(name, prefix=None, suffix=None):
    """Remove `prefix` and / or `suffix` from `name`.

    Examples
    --------
    >>> strip_ps('good_god_gomer', 'good', 'gomer')
    '_god_'
    """
    if prefix and name.startswith(prefix):
        name = name[len(prefix):]
    if suffix and name.endswith(suffix):
        name = name[:-len(suffix)]
    return name


def format_unprintable(string):
    r"""Escape tabs (\t), newlines (\n), etc. for system commands and printing.

    Examples
    --------
    >>> format_unprintable('\t')
    '\\t'
    """
    return repr(string).strip("'")


# =============================================================================
# Filesystem operation
# =============================================================================

@contextmanager
def switch_paths(working_path):
    current_path = os.getcwd()
    try:
        os.chdir(working_path)
        yield
    except:
        raise
    finally:
        os.chdir(current_path)


def copyfile(infile, outfile, mode=None):
    shutil.copyfile(infile, outfile)
    if mode is not None:
        os.chmod(outfile, mode)


def makedirs(path, mode=None, exist_ok=True):
    if mode is None:
        os.makedirs(path, exist_ok=exist_ok)
    else:
        # Don't think this works as expected...
        original_umask = os.umask(0)
        try:
            os.makedirs(path, mode=mode, exist_ok=exist_ok)
        finally:
            os.umask(original_umask)


def make_tarfile(source_dir, output_filename):
    """Compress folder into a `*.tar.gz` file."""
    import tarfile
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=op.basename(source_dir))
