import io
import logging
import shlex
import subprocess
from pathlib import Path
from typing import Callable, Union

import kmbio.PDB
import modeller
from Bio.AlignIO import MultipleSeqAlignment
from modeller.automodel import assess, automodel, autosched

from kmtools import py_tools, system_tools

logger = logging.getLogger(__name__)


def run_modeller(structure, alignment, temp_dir: Union[str, Path, Callable]):
    """Run Modeller to create a homology model.

    Args:
        structure: Structure of the template protein.
        alignment_file: Alignment of the target sequence(s) to chain(s) of the template structure.
        temp_dir: Location to use for storing Modeller temporary files and output.

    Returns:
        results: List of successfully calculated models.
    """
    if isinstance(structure, (str, Path)):
        structure = kmbio.PDB.load(structure)

    if callable(temp_dir):
        temp_dir = Path(temp_dir())
    else:
        temp_dir = Path(temp_dir)

    assert len(alignment) == 2
    target_id = alignment[0].id
    template_id = alignment[1].id

    kmbio.PDB.save(structure, temp_dir.joinpath(f"{template_id}.pdb"))
    alignment_file = temp_dir.joinpath(f"{template_id}-{target_id}.aln")
    write_pir_alignment(alignment, alignment_file)

    # Don't display log messages
    modeller.log.none()

    # Create a new MODELLER environment
    env = modeller.environ()

    # Directories for input atom files
    env.io.atom_files_directory = [str(temp_dir)]
    env.schedule_scale = modeller.physical.values(default=1.0, soft_sphere=0.7)

    # Selected atoms do not feel the neighborhood
    # env.edat.nonbonded_sel_atoms = 2
    env.io.hetatm = True  # read in HETATM records from template PDBs
    env.io.water = True  # read in WATER records (including waters marked as HETATMs)

    a = automodel(
        env,
        # alignment filename
        alnfile=str(alignment_file),
        # codes of the templates
        knowns=(str(template_id)),
        # code of the target
        sequence=str(target_id),
        # wich method for validation should be calculated
        assess_methods=(assess.DOPE, assess.normalized_dope, assess.GA341),
    )
    a.starting_model = 1  # index of the first model
    a.ending_model = 1  # index of the last model

    # Very thorough VTFM optimization:
    a.library_schedule = autosched.slow
    a.max_var_iterations = 300

    # Thorough MD optimization:
    # a.md_level = refine.slow
    a.md_level = None

    # Repeat the whole cycle 2 times and do not stop unless obj.func. > 1E6
    # a.repeat_optimization = 2

    a.max_molpdf = 2e5

    with py_tools.log_print_statements(logger), system_tools.switch_paths(temp_dir):
        a.make()

    assert len(a.outputs) == 1
    return a.outputs[0]


def write_pir_alignment(alignment: MultipleSeqAlignment, file: Path):
    """Write `alignment` into a pir file."""
    assert len(alignment) == 2
    seqrec_1, seqrec_2 = alignment
    with open(file, "wt") as fout:
        # Sequence
        fout.write(f">P1;{seqrec_1.id}\n")
        fout.write(f"sequence:{seqrec_1.id}:.:.:.:.::::\n")
        fout.write(f"{seqrec_1.seq}*\n\n")
        # Structure
        fout.write(f">P1;{seqrec_2.id}\n")
        fout.write(f"structure:{seqrec_2.id}:.:.:.:.::::\n")
        fout.write(f"{seqrec_2.seq}*\n")


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


def protonate(method, structure, **kwargs):
    ...
