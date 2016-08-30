import os.path as op
import pandas as pd
import kmtools
import pytest


@pytest.mark.parametrize("uniprot_id, uniprot_sequence", [
    ('P42212', """\
MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL
VTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLV
NRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLAD
HYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK
""".replace('\n', '').replace(' ', '')),
])
def test_get_uniprot_sequence(uniprot_id, uniprot_sequence):
    assert str(kmtools.sequence_tools.get_uniprot_sequence(uniprot_id).seq) == uniprot_sequence


def test_get_crossmapping():
    sequence_pairs_df = pd.read_csv(
        op.join(op.splitext(__file__)[0], 'test_alignment.tsv.gz'),
        sep='\t')

    sequence_pairs_df['alignment_ref'], sequence_pairs_df['alignment_alt'] = (
        list(zip(*sequence_pairs_df.apply(
            lambda x: kmtools.sequence_tools.align_pairwise(
                x['sequence_ref'], x['sequence_alt']), axis=1))))
    assert (sequence_pairs_df['alignment_ref'].str.replace('-', '').str.len() ==
            sequence_pairs_df['sequence_ref'].str.len()).all().all()
    assert (sequence_pairs_df['alignment_alt'].str.replace('-', '').str.len() ==
            sequence_pairs_df['sequence_alt'].str.len()).all().all()
    assert (sequence_pairs_df['alignment_ref'].str.len() ==
            sequence_pairs_df['alignment_alt'].str.len()).all().all()

    sequence_pairs_df['mapping_ref'], sequence_pairs_df['mapping_alt'] = (
        list(zip(*sequence_pairs_df.apply(
            lambda x: kmtools.sequence_tools.get_crossmapping(
                x['alignment_ref'], x['alignment_alt']), axis=1))))
    assert ((sequence_pairs_df['mapping_ref'].str.count(',') + 1) ==
            (sequence_pairs_df['sequence_alt'].str.len())).all().all()
    assert ((sequence_pairs_df['mapping_alt'].str.count(',') + 1) ==
            (sequence_pairs_df['sequence_ref'].str.len())).all().all()

    def validate_mapping(sequence_ref, sequence_alt, mapping_ref2alt):
        return all(
            (sequence_ref[int(x) - 1] == sequence_alt[i])
            for (i, x)
            in enumerate(mapping_ref2alt.split(','))
            if x
        )
    assert sequence_pairs_df[['sequence_ref', 'sequence_alt', 'mapping_ref']].apply(
        lambda x: validate_mapping(*x), axis=1
    ).all()
    assert sequence_pairs_df[['sequence_alt', 'sequence_ref', 'mapping_alt']].apply(
        lambda x: validate_mapping(*x), axis=1
    ).all()
