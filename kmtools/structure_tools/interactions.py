import logging

import pandas as pd

from kmbio.PDB import NeighborSearch
from kmtools import df_tools
from kmtools.structure_tools import AAA_DICT, COMMON_HETATMS

logger = logging.getLogger(__name__)


def _iter_interchain_ns(structure, interchain=True):
    for model_1_idx, model_1 in enumerate(structure):
        for chain_1_idx, chain_1 in enumerate(model_1):
            atom_list = list(chain_1.atoms)
            if not atom_list:
                logger.debug("Skipping %s %s because it is empty.", model_1, chain_1)
                continue
            chain_1_ns = NeighborSearch(atom_list)
            if not interchain:
                yield (model_1_idx, model_1, chain_1_idx, chain_1, chain_1_ns, model_1_idx,
                       model_1, chain_1_idx, chain_1)
            else:
                for model_2_idx, model_2 in enumerate(structure):
                    if model_1_idx > model_2_idx:
                        continue
                    for chain_2_idx, chain_2 in enumerate(model_2):
                        if model_1_idx == model_2_idx and chain_1_idx > chain_2_idx:
                            continue
                        yield (model_1_idx, model_1, chain_1_idx, chain_1, chain_1_ns, model_2_idx,
                               model_2, chain_2_idx, chain_2)


def _get_interactions(structure, r_cutoff, interchain):
    results = []
    for (model_1_idx, model_1, chain_1_idx, chain_1, chain_1_ns, model_2_idx, model_2, chain_2_idx,
         chain_2) in _iter_interchain_ns(
             structure, interchain=interchain):
        chain_1_residue_ids = [r.id for r in chain_1]
        for residue_2_idx, residue_2 in enumerate(chain_2):
            if residue_2.resname in COMMON_HETATMS:
                continue
            seen = set()
            chain_1_interacting_residues = [
                r for a in residue_2 for r in chain_1_ns.search(a.coord, r_cutoff, 'R')
                if r.id not in seen and not seen.add(r.id)
            ]
            for residue_1 in chain_1_interacting_residues:
                if residue_1.resname in COMMON_HETATMS:
                    continue
                residue_1_idx = chain_1_residue_ids.index(residue_1.id)
                row = (structure.id, model_1.id, model_2.id, chain_1.id, chain_2.id, residue_1_idx,
                       residue_2_idx, residue_1.id[1], residue_2.id[1], residue_1.resname,
                       residue_2.resname, AAA_DICT.get(residue_1.resname),
                       AAA_DICT.get(residue_2.resname))
                results.append(row)


def _gen_interaction_df(interactions):
    columns = [
        'structure_id',
        'model_id_1',
        'model_id_2',
        'chain_id_1',
        'chain_id_2',
        'residue_idx_1',
        'residue_idx_2',
        'residue_id_1',
        'residue_id_2',
        'residue_name_1',
        'residue_name_2',
        'residue_aa_1',
        'residue_aa_2',
    ]
    df = pd.DataFrame(interactions, columns=columns)
    df_reverse = df.rename(columns=df_tools.get_reverse_column).reindex(columns=pd.Index(columns))
    _length_before = len(df)
    df = pd.concat([df, df_reverse]).drop_duplicates()
    assert len(df) == _length_before * 2  # Make sure that there were no duplicates
    df = df.sort_values(columns)  # Will not work well for uppercase / lowercase columns
    df.index = range(len(df))
    return df


def get_interactions(structure, r_cutoff=5.0, interchain=True):
    """Return residue-residue interactions within and between chains in `structure`.

    .. todo::

        This could probably be sped up by using the
        :py:meth:`kmbio.PDB.NeighborSearch.search_all` method.

    Parameters
    ----------
    structure : kmbio.PDB.Structure.Structure
        Structure to analyse.
    interchain : :class:`bool`
        Whether to include interactions between chains.

    Notes
    -----
    - For each chain, `residue.id[1]` is unique.
    """
    interactions = _get_interactions(structure, r_cutoff, interchain)
    interactions_df = _gen_interaction_df(interactions)
    return interactions_df
