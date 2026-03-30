import io

import pytest

from string_gsea.get_species import GetTaxonID, OxFieldsZip  # Import GetTaxonID
from string_gsea.gsea_utilities import get_rank_files  # Import get_rank_files

# Expected species ID (adjust if necessary based on fasta_test.zip content)
EXPECTED_SPECIES_ID = 559292  # This was for the fasta test
# Expected OX values (adjust if necessary) - assuming one fasta file inside with one OX=9606 entry
EXPECTED_STRING_SPECIES_ID = 4932


def test_get_species_from_oxes(fasta_test_zip):
    """Tests the get_species_from_oxes function."""
    if not fasta_test_zip.exists():
        pytest.fail(f"Test zip file not found at {fasta_test_zip}")

    species_id = OxFieldsZip.get_species_from_oxes(fasta_test_zip)
    assert species_id == EXPECTED_SPECIES_ID, (
        f"Expected species ID {EXPECTED_SPECIES_ID}, but got {species_id}"
    )


def test_get_ox_fields():
    fasta_content = (
        b">sp|P0DTC2|SPIKE_SARS2 Spike glycoprotein OS=Severe acute respiratory syndrome coronavirus 2 OX=2697049 GN=S\n"
        b">sp|P0DTC1|NCAP_SARS2 Nucleoprotein OS=Severe acute respiratory syndrome coronavirus 2 OX=2697049 GN=N\n"
    )
    expected_ox = [2697049, 2697049]
    fasta_stream = io.BytesIO(fasta_content)
    result = OxFieldsZip.get_ox_fields(fasta_stream)
    assert result == expected_ox


@pytest.mark.integration
def test_GetTaxonID_determine_species(yeast_rnk_zip) -> None:
    """Tests GetTaxonID.get_species_from_rank_file using a yeast rank file."""
    if not yeast_rnk_zip.exists():
        pytest.fail(f"Test rank zip file not found at {yeast_rnk_zip}")

    # Use get_rank_files to read the data
    rank_dataframes = get_rank_files(yeast_rnk_zip)

    # Check if any rank files were found and loaded
    if not rank_dataframes:
        pytest.fail(f"No .rnk files found or loaded from {yeast_rnk_zip}")

    # Get the first DataFrame from the dictionary (assuming at least one exists)
    # The key doesn't matter for this test, only the content of the DataFrame
    first_df = next(iter(rank_dataframes.values()))

    # Call the function under test
    # Note: This test makes live calls to the STRING API
    try:
        species_id = GetTaxonID.determine_species(
            first_df, nr=10
        )  # Use fewer samples for faster test
        assert species_id == EXPECTED_STRING_SPECIES_ID, (
            f"Expected species ID {EXPECTED_STRING_SPECIES_ID}, but got {species_id}"
        )
    except ValueError as e:
        pytest.fail(
            f"GetTaxonID.determine_species raised an unexpected ValueError: {e}"
        )
    except Exception as e:
        pytest.fail(f"An unexpected error occurred during API call or processing: {e}")
