"""Legacy code to support the ELASPIC pipeline until we move to ELASPIC v2."""
from kmtools.structure_tools import AAA_DICT


def get_chain_sequence_and_numbering(chain, domain_def_tuple=None, include_hetatms=False):
    """Get the amino acid sequence and a list of residue ids for the given chain.

    Parameters
    ----------
    chain : Bio.PDB.Chain.Chain
        The chain for which to get the amino acid sequence and numbering.
    """
    if domain_def_tuple is not None:
        start_resid, end_resid = domain_def_tuple

    chain_numbering = []
    chain_numbering_extended = []
    chain_sequence = []
    inside_domain = False
    for res in chain:
        #
        resid = str(res.id[1]) + res.id[2].strip()

        if domain_def_tuple is None or resid == start_resid:
            inside_domain = True

        if inside_domain and (include_hetatms or res.resname in AAA_DICT):
            chain_numbering.append(res.id[1])
            chain_numbering_extended.append(resid)
            chain_sequence.append(AAA_DICT.get(res.resname, '.'))

        if domain_def_tuple is not None and resid == end_resid:
            inside_domain = False

    chain_sequence = ''.join(chain_sequence)
    return chain_sequence, chain_numbering_extended
