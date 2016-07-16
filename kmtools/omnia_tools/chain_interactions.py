"""Find all pairs of interacting amino acids using `mdtraj`.

.. note::

    Work in progress...

"""
import os.path as op
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display
import mdtraj as md
# import kmtools.pdb_tools.sifts

PDB_CACHE_DIR = None


def main(pdb_file):
    t = md.load(op.join(PDB_CACHE_DIR, pdb_file))
    distances, residue_pairs = md.compute_contacts(t, [[0, i] for i in range(700)])
    sum(r.is_protein for r in t.topology.residues)
    md.compute_contacts(t, [[0, ]])
    t.topology.chain(0).topology.to_dataframe()[0]
    t.topology.to_fasta()

    r = t.topology.residue(0)
    r.resSeq
    r

    for r in t.topology.residues:
        print(r.chain.index, r.index, r.resSeq, r.is_protein)

    distances, residue_pairs = md.compute_contacts(t)
    residue_pairs.shape
    plt.hist(distances[0])
    assert distances.shape[1] == residue_pairs.shape[0]
    t.topology.residue(770)

    contacts_df = pd.concat([
        pd.Series(residue_pairs[:, 0], name='residue_1_idx'),
        pd.Series(residue_pairs[:, 1], name='residue_1_idx'),
        pd.Series(distances[0], name='distance')
    ], axis=1)
    contacts_df.head()

    a, b = t.topology.to_dataframe()
    display(a.tail())
    b

    # g = t.topology.to_bondgraph()

    t.topology.chain(2).topology.to_dataframe()[0]

    contacts_df.head()

    contacts_df.max()
    # sifts_df = kmtools.pdb_tools.sifts.get_sifts_data('1CSE')


if __name__ == '__main__':
    main('1cse.pdb')
