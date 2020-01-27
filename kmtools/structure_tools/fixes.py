import contextlib
import logging
import shlex
import subprocess
from pathlib import Path
from typing import IO, Union

logger = logging.getLogger(__name__)


def protonate(
    input_file: Union[str, Path, IO], output_file: Union[str, Path, IO], method: str = "openmm"
) -> None:
    """Add hydrogen atoms to PDB structure.

    Args:
        input_file: Input PDB file or IO object.
        output_file: Output PDB file or IO object.
        method: Method to use for adding hydrogens. Supported methods are "openmm" and "reduce".
    """
    supported_methods = ["openmm", "reduce", "reduce-FLIP", "reduce-NOFLIP", "reduce-BUILD"]
    if method not in supported_methods:
        raise ValueError(f"The only supported methods are '{supported_methods}'.")
    with contextlib.ExitStack() as stack:
        if isinstance(input_file, (str, Path)):
            input_file = stack.enter_context(open(input_file, "rt"))
        if isinstance(output_file, (str, Path)):
            output_file = stack.enter_context(open(output_file, "wt"))
        if method.startswith("openmm"):
            _protonate_with_openmm(input_file, output_file)
        elif method.startswith("reduce"):
            _protonate_with_reduce(input_file, output_file, method="")


def _protonate_with_reduce(input_file: IO, output_file: IO, method: str = "") -> None:
    assert method in ["", "FLIP", "NOFLIP", "BUILD"]
    system_command = "reduce {} -".format("-" + method.upper() if method else "")
    proc = subprocess.run(
        shlex.split(system_command),
        input=input_file.read(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    if proc.returncode != 1:
        logger.warning(
            f"Command '{system_command}' returned non-zero exit status {proc.returncode}."
        )
    output_file.write(proc.stdout)


def _protonate_with_openmm(input_file: IO, output_file: IO) -> None:
    from simtk.openmm import app

    pdb = app.PDBFile(input_file)
    modeller = app.Modeller(pdb.getTopology(), pdb.getPositions())
    modeller.addHydrogens(pH=7.0)
    app.PDBFile.writeFile(modeller.topology, modeller.positions, output_file, keepIds=True)
