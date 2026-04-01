"""Tests for species detection — FASTA OX parsing and STRING API lookup."""

import io

import pytest

from string_gsea.gsea_utilities import get_rank_files
from string_gsea.species_detection import get_ox_fields, get_species_from_oxes
from string_gsea.string_api_client import determine_species_from_identifiers

# Expected species ID (adjust if necessary based on fasta_test.zip content)
EXPECTED_SPECIES_ID = 559292  # This was for the fasta test
# Expected OX values (adjust if necessary) - assuming one fasta file inside with one OX=9606 entry
EXPECTED_STRING_SPECIES_ID = 4932


def test_get_species_from_oxes(fasta_test_zip):
    """Tests the get_species_from_oxes function."""
    if not fasta_test_zip.exists():
        pytest.fail(f"Test zip file not found at {fasta_test_zip}")

    species_id = get_species_from_oxes(fasta_test_zip)
    assert species_id == EXPECTED_SPECIES_ID, f"Expected species ID {EXPECTED_SPECIES_ID}, but got {species_id}"


def test_get_ox_fields():
    fasta_content = (
        b">sp|P0DTC2|SPIKE_SARS2 Spike glycoprotein"
        b" OS=Severe acute respiratory syndrome coronavirus 2 OX=2697049 GN=S\n"
        b">sp|P0DTC1|NCAP_SARS2 Nucleoprotein"
        b" OS=Severe acute respiratory syndrome coronavirus 2 OX=2697049 GN=N\n"
    )
    expected_ox = [2697049, 2697049]
    fasta_stream = io.BytesIO(fasta_content)
    result = get_ox_fields(fasta_stream)
    assert result == expected_ox


@pytest.mark.integration
def test_determine_species(yeast_rnk_zip) -> None:
    """Tests determine_species_from_identifiers using a yeast rank file."""
    if not yeast_rnk_zip.exists():
        pytest.fail(f"Test rank zip file not found at {yeast_rnk_zip}")

    rank_lists = get_rank_files(yeast_rnk_zip)

    if not rank_lists:
        pytest.fail(f"No .rnk files found or loaded from {yeast_rnk_zip}")

    first_rl = next(iter(rank_lists.values()))
    identifiers = list(first_rl.entries.keys())

    try:
        species_id = determine_species_from_identifiers(identifiers, nr=10)
        assert species_id == EXPECTED_STRING_SPECIES_ID, (
            f"Expected species ID {EXPECTED_STRING_SPECIES_ID}, but got {species_id}"
        )
    except ValueError as e:
        pytest.fail(f"determine_species_from_identifiers raised an unexpected ValueError: {e}")
    except Exception as e:
        pytest.fail(f"An unexpected error occurred during API call or processing: {e}")
