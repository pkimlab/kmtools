import elaspic.helper
from . import pdb_tools


def launch_pymol(quiet=False, gui=False):
    """Launch PyMOL.

    This code has to be inside a function in order not to break doctest.
    """
    import __main__
    args = ''
    if quiet:
        args += 'q'
    if gui:
        args += 'c'
    if args:
        __main__.pymol_argv = ['pymol', '-' + args]
    import pymol
    pymol.finish_launching()
    return pymol


def show_pymol(domain_definition, list_of_mutations=None):
    """Colour residues by their sasa score instead of the b factor."""
    pymol = launch_pymol()
    d, t, m = domain_definition
    path_to_pdb = d.path_to_data + m.model_filename
    structure_name = m.model_filename.strip('.pdb')

    sasa_score = [float(score) for score in m.sasa_score.split(',')]

    # Load the pdb
    pymol.cmd.load(path_to_pdb, structure_name)

    # clear out the old B Factors
    pymol.cmd.alter('%s' % structure_name, 'b=0.0')

    # update the B Factors with new properties
    for idx, sasa in enumerate(sasa_score):
        resnum = idx + 1
        pymol.cmd.alter('%s and resi %i' % (structure_name, resnum), 'b=%f' % sasa)

    #    color the protein based on the new B Factors
    pymol.cmd.spectrum("b", "protA and n. CA")

    if list_of_mutations:
        pymol.stored.structure_aaa = None
        for pymol.stored.mutation in list_of_mutations:

            from_aa = pymol.stored.mutation[0]
            # to_aa = pymol.stored.mutation[-1]
            pymol.stored.uniprot_position = int(pymol.stored.mutation[1:-1])
            resnum = (
                pymol.stored.uniprot_position -
                elaspic.helper.decode_domain_def(t.domain_def)[0] + 1
            )

            pymol.stored.mutation_a = from_aa.upper()
            pymol.cmd.iterate(
                '%s and resi %i and n. CA' % (structure_name, resnum),
                'print resn; stored.structure_aaa = resn')
            assert pdb_tools.AAA_DICT[pymol.stored.structure_aaa] == pymol.stored.mutation_a
            pymol.cmd.label(
                '%s and resi %i and n. CA' % (structure_name, resnum), "stored.mutation")

    pymol.cmd.hide('all')
    pymol.cmd.show('cartoon')
    pymol.cmd.show('labels')

    pymol.cmd.set('surface_quality', 1)
    pymol.cmd.set('transparency', 0.6)

    pymol.cmd.refresh()
#    pymol.cmd.save(path_to_pdb.replace('.pdb','-sasa_bfactors.pdb'))
