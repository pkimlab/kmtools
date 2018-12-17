import io
import logging
import shlex
import subprocess
from pathlib import Path
from textwrap import dedent
from typing import Mapping

import Bio.SeqIO

logger = logging.getLogger(__name__)


def run_psipred(fasta_file: Path, hhblits_database: Path, nproc: int = None) -> Path:
    """Predict secondary structure of sequence in `fasta_file`.

    Args:
        fasta_file: File contianing the input sequence.
        hhblits_database: Root folder of the ``hhblits`` database
            (e.g. ``./uniprot20_2016_02``, not ``./uniprot20_2016_02/uniprot20_2016_02``).
        nproc: Number of cores that should be used.

    Returns:
        Path of the file containing predicted secondary structure.
    """
    if nproc is not None:
        env = {"OMP_NUM_THREADS": str(nproc)}
    else:
        env = {}

    hhblits_system_command = dedent(
        f"""\
    hhblits \
        -i {fasta_file} \
        -oa3m {fasta_file.with_suffix('.a3m')} \
        -all \
        -maxfilt 100000 \
        -realign_max 100000 \
        -B 100000 \
        -Z 100000 \
        -d {hhblits_database}/{hhblits_database.name}
    """
    )
    _execute(hhblits_system_command, env)

    hhfilter_system_command = dedent(
        f"""\
    hhfilter \
        -i {fasta_file.with_suffix('.a3m')} \
        -o {fasta_file.with_suffix('.filt.a3m')} \
        -id 90 \
        -neff 15 \
        -qsc -30
    """
    )
    _execute(hhfilter_system_command, env)

    addss_system_command = dedent(
        f"""\
    addss.pl \
        {fasta_file.with_suffix('.filt.a3m')} \
        {fasta_file.with_suffix('.filt.withss.a3m')} \
        -a3m
    """
    )
    _execute(addss_system_command, env)

    return fasta_file.with_suffix(".filt.withss.a3m")


def read_psipred(psipred_file: Path) -> str:
    """Read fasta file produced by `run_psipred`."""
    buf = io.StringIO()
    with psipred_file.open("rt") as fin:
        for i, line in enumerate(fin):
            buf.write(line)
            if i >= 1:
                break
    buf.seek(0)
    psipred = Bio.SeqIO.read(buf, format="fasta")
    return str(psipred.seq)


def _execute(system_command: str, env: Mapping[str, str] = None) -> None:
    env = env if env is not None else {}
    logger.debug(f"Running system command: '{system_command}'...")
    proc = subprocess.run(
        shlex.split(system_command),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        env=env,
    )
    if proc.returncode:
        logger.error(proc.stdout)
        logger.error(proc.stderr)
        raise subprocess.CalledProcessError(
            proc.returncode, system_command, proc.stdout, proc.stderr
        )
