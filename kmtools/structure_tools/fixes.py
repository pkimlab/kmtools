import io
import shlex
import subprocess

import kmbio.PDB


def protonate(method, structure, **kwargs):
    ...


def _protonate_with_reduce(structure, method=None):
    assert method in [None, "FLIP", "NOFLIP", "BUILD"]
    system_command = "reduce {} -".format("-" + method if method is not None else "")
    input_file = io.StringIO()
    kmbio.PDB.save(structure, input_file)
    input_file.seek(0)
    output_file = io.StringIO()
    proc = subprocess.run(
        shlex.split(system_command),
        input=input_file,
        stdout=output_file,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    if proc.returncode not in [0, 1]:
        raise subprocess.CalledProcessError(proc)


def _protonate_with_openmm(structure, num_steps=None):
    ...
