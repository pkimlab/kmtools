import logging
import shlex
import subprocess
from contextlib import closing
from pathlib import Path

from kmtools import py_tools

from .types import AddSS, HHBlits, HHFilter, HHInput, HHMake, HHMakeModel, HHSearch

logger = logging.getLogger(__name__)


def run(cmd: str, cwd: Path) -> None:
    logger.debug("system_commmand: '%s'", cmd)
    with closing(py_tools.LogPipe(logger.debug)) as pipe:
        subprocess.run(  # type: ignore
            shlex.split(cmd), stdout=pipe, stderr=pipe, universal_newlines=True, cwd=cwd, check=True
        )


def hhblits(input: HHInput, database: Path, extra_args: str = "") -> HHBlits:
    """Run ``hhblits``, HMM-HMM-based lightning-fast iterative sequence search.

    Args:
        input: :any:`HHInput` object containing required input data.
        database: Database used for constructing the ``hhblits`` alignment.
        extra_args: Extra parameters to pass down to the executable.

            - For homology modeling (hhpred), use:
              ``-n 3 -mact 0.5``.
            - For secondary structure / coevolution prediction, use:
              ``-all -maxfilt 100000 -realign_max 100000 -B 100000 -Z 100000``.

    Returns:
        :any:`HHBlits` object containing results.
    """
    output_file_base = input.temp_dir.joinpath(input.sequence_file.name)
    data = HHBlits(
        a3m_file=output_file_base.with_suffix(".a3m"),
        hhm_file=output_file_base.with_suffix(".hhm"),
        hhblits_hhr_file=output_file_base.with_suffix(".hhr"),
        hhblits_database=database,
        hhblits_extra_args=extra_args,
        **vars(input),
    )
    cmd = (
        f"hhblits"
        f" -i '{data.sequence_file}'"
        f" -o '{data.hhblits_hhr_file}'"
        f" -oa3m '{data.a3m_file}'"
        f" -ohhm '{data.hhm_file}'"
        f" -oalis '{data.temp_dir.joinpath('hhblits_alis.a3m')}'"
        f" -d '{data.hhblits_database}/{data.hhblits_database.name}'"
        f" -M first"
        f" {data.hhblits_extra_args}"
    )
    run(cmd, data.temp_dir)
    return data


def hhfilter(input: HHBlits, extra_args: str = "") -> HHFilter:
    """Filter an alignment using ``hhfilter``.

    Args:
        input: :any:`HHBlits` object containing required input data.
        extra_args: Extra parameters to pass down to the executable.

            - For secondary structure / coevolution prediction prediction, use:
              ``-id 90 -neff 15 -qsc -30``.

    Returns:
        :any:`HHFilter` object containing results.
    """
    data = HHFilter(hhfilter_extra_args=extra_args, **vars(input))
    data.a3m_file = input.a3m_file.with_suffix(".filt.a3m")
    cmd = (
        f"hhfilter"
        f" -i '{input.a3m_file}'"
        f" -o '{data.a3m_file}'"
        f" {data.hhfilter_extra_args}"
    )
    run(cmd, data.temp_dir)
    return data


def addss(input: HHBlits, extra_args: str = "") -> AddSS:
    """Add secondary structure information to an alignment using ``addss.pl``.

    Args:
        input: :any:`HHBlits` object containing required input data.
        extra_args: Extra parameters to pass down to the executable.

    Returns:
        :any:`AddSS` object containing results.
    """
    data = AddSS(addss_extra_args=extra_args, **vars(input))
    data.a3m_file = input.a3m_file.with_suffix(".withss.a3m")
    cmd = (
        #
        f"addss.pl"
        f" '{input.a3m_file}'"
        f" '{data.a3m_file}'"
        f" -a3m"
        f" {data.addss_extra_args}"
    )
    run(cmd, data.temp_dir)
    return data


def hhmake(input: HHBlits, extra_args: str = "") -> HHMake:
    """Run ``hhmake`` to build an HMM from an input alignment.

    Args:
        input: :any:`HHBlits` object containing required input data.
        extra_args: Extra parameters to pass down to the executable.

    Returns:
        :any:`HHMake` object containing results.
    """
    data = HHMake(hhmake_extra_args=extra_args, **vars(input))
    data.hhm_file = input.a3m_file.with_suffix(".hhm")
    cmd = (
        #
        f"hhmake"
        f" -i '{data.a3m_file}'"
        f" -o '{data.hhm_file}'"
        f" {data.hhmake_extra_args}"
    )
    run(cmd, data.temp_dir)
    return data


def hhsearch(input: HHMake, database: Path, extra_args: str = "") -> HHSearch:
    """Run ``hhsearch`` to "search a database of HMMs with a query alignment or query HMM".

    Args:
        input: :any:`HHMake` object containing required input data.
        database: Database used for constructing the ``hhblits`` alignment.
        extra_args: Extra parameters to pass down to the executable.

            - For homology modelling, use: ``-mact 0.05 -e 0.1 -glob``.

    Returns:
        :any:`HHSearch` object containing results.
    """
    data = HHSearch(
        hhsearch_hhr_file=input.hhm_file.with_suffix(".hhsearch.hhr"),
        hhsearch_tab_file=input.hhm_file.with_suffix(".hhsearch.tab"),
        hhsearch_database=database,
        hhsearch_extra_args=extra_args,
        **vars(input),
    )
    cmd = (
        f"hhsearch"
        f" -i '{data.hhm_file}'"
        f" -o '{data.hhsearch_hhr_file}'"
        f" -d '{data.hhsearch_database}/{data.hhsearch_database.name}'"
        f" -atab '{data.hhsearch_tab_file}'"
        f" {data.hhsearch_extra_args}"
    )
    run(cmd, data.temp_dir)
    return data


def hhmakemodel(input: HHSearch, extra_args: str = "-v 1 -m 1") -> HHMakeModel:
    """Run ``hhmakemodel.pl``.

    Args:
        input: :any:`HHSearch` object containing required input data.
        extra_args: Extra parameters to pass down to the executable.

    Returns:
        :any:`HHMakeModel` object containing results.
    """
    data = HHMakeModel(
        hhmakemodel_pir_file=input.temp_dir.joinpath("hhmakemodel_pir_file"),
        hhmakemodel_extra_args=extra_args,
        **vars(input),
    )
    cmd = (
        f"hhmakemodel.pl"
        f" '{data.hhsearch_hhr_file}'"
        f" -q '{data.a3m_file}'"
        f" -d '{data.hhsearch_database}/{data.hhsearch_database.name}'"
        f" -pir '{data.hhmakemodel_pir_file}'"
        f" {data.hhmakemodel_extra_args}"
    )
    run(cmd, data.temp_dir)
    return data
