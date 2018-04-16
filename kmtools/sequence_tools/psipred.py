import io
import logging
import shlex
import subprocess
from pathlib import Path
from textwrap import dedent

import Bio.SeqIO

logger = logging.getLogger(__name__)


def run_psipred(fasta_file: Path) -> Path:
    """Predict secondary structure of sequence in `fasta_file`."""
    hhblits_system_command = dedent(f"""\
    hhblits \\
        -i {fasta_file} \\
        -oa3m {fasta_file.with_suffix('.a3m')} \\
        -all \\
        -maxfilt 100000 \\
        -realign_max 100000 \\
        -B 100000 \\
        -Z 100000 \\
        -d /home/kimlab1/database_data/hh-suite/uniprot20_2016_02/uniprot20_2016_02
    """)
    _execute(hhblits_system_command)

    hhfilter_system_command = dedent(f"""\
    hhfilter \\
        -i {fasta_file.with_suffix('.a3m')} \\
        -o {fasta_file.with_suffix('.filt.a3m')} \\
        -id 90 \\
        -neff 15 \\
        -qsc -30
    """)
    _execute(hhfilter_system_command)

    addss_system_command = dedent(f"""\
    addss.pl \\
        {fasta_file.with_suffix('.filt.a3m')} \\
        {fasta_file.with_suffix('.filt.withss.a3m')} \\
        -a3m
    """)
    _execute(addss_system_command)

    return fasta_file.with_suffix('.filt.withss.a3m')


def read_psipred(fasta_withss_file: Path) -> Path:
    """Read fasta file produced by `run_psipred`."""
    buf = io.StringIO()
    with fasta_withss_file.open('rt') as fin:
        for i, line in enumerate(fin):
            buf.write(line)
            if i >= 1:
                break
    buf.seek(0)
    ss = Bio.SeqIO.read(buf, format='fasta')
    return ss


def _execute(system_command: str) -> None:
    logger.debug(f"Running system command: '{system_command}'...")
    proc = subprocess.run(
        shlex.split(system_command),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        check=True,
    )
    if proc.returncode:
        raise subprocess.CalledProcessError(proc.returncode, system_command, proc.stdout,
                                            proc.stderr)
