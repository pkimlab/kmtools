from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile


# %%
def fix_pdb(input_file, output_file):
    """
    .. note::
        Set ``export OPENMM_CPU_THREADS=1`` if you want this function
        to use only one thead.
    """
    fixer = PDBFixer(filename=input_file)
    fixer.findMissingResidues()
    # Replace non-standard residues with their standard equiv.
    # (But if atoms have to be added, they get added using the `.addMissingAtoms()` command)
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    # Find missing heavy atoms
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    # fixer.removeHeterogens(True)
    fixer.addMissingHydrogens(7.0)
    with open(output_file, 'w') as ofh:
        PDBFile.writeFile(fixer.topology, fixer.positions, ofh, keepIds=True)
