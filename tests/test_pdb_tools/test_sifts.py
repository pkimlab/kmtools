import pytest
import kmtools


@pytest.mark.parametrize("pdb_id,pdb_chains,pdb_mutations", [
    ('2qja', 'C', 'T74A,L78M,G79K,L80Y'),
])
def test_sifts_exception(pdb_id, pdb_chains, pdb_mutations):
    """Test the case where PDB AA and UniProt AA are different."""
    sifts_df = kmtools.pdb_tools.sifts.get_sifts_data(pdb_id)
    with pytest.raises(kmtools.pdb_tools.sifts.SIFTSError):
        kmtools.pdb_tools.sifts.convert_pdb_mutations_to_uniprot(
            pdb_id, pdb_chains, pdb_mutations, sifts_df=sifts_df)
