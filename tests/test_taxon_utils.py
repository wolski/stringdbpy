import pytest
import polars as pl
from pathlib import Path

# Import the TaxonUtils class from your module
from string_gsea.get_species import TaxonUtils

# Create a fixture for the TaxonUtils instance
@pytest.fixture
def taxon_utils():
    return TaxonUtils()

# Check file existence for skipping tests if necessary
species_file_exists = Path(__file__).resolve().parent.parent / 'data' / 'mappings' / "species.v12.0.zip"
ncbi_file_exists = Path(__file__).resolve().parent.parent / 'data' / 'mappings' / "NCBI_nodes.zip"

@pytest.mark.skipif(not species_file_exists, reason="Species zip file not found")
def test_read_species_string_data(taxon_utils):
    # Access the private method using the instance
    species_df = taxon_utils._read_species_string_data()
    assert species_df is not None, "The species DataFrame should not be None."

    # Check that the '#taxon_id' column was renamed to "taxon_id"
    assert "taxon_id" in species_df.columns, "Expected column 'taxon_id' was not found."
    assert "#taxon_id" not in species_df.columns, "Column '#taxon_id' was not renamed."

    # Optionally, ensure that some rows are loaded.
    assert species_df.height > 0, "No rows found in species data."


@pytest.mark.skipif(not ncbi_file_exists, reason="NCBI zip file not found")
def test_read_ncbi_nodes_data(taxon_utils):
    # Access the private method using the instance
    ncbi_df = taxon_utils._read_ncbi_nodes_data()
    assert ncbi_df is not None, "The NCBI DataFrame should not be None."

    # Check that both required columns are present.
    for col in ("taxon_id", "parent_taxon_id"):
        assert col in ncbi_df.columns, f"Expected column '{col}' not found in NCBI data."

    # Optionally, check that the file has some rows.
    assert ncbi_df.height > 0, "No rows found in NCBI data."


@pytest.mark.skipif(not (species_file_exists and ncbi_file_exists), reason="Necessary data files not found")
def test_get_organism_for_string(taxon_utils):
    # --- Test case 1: Direct match ---
    # Here, we assume that taxon id 559292 is present in species data.
    result = taxon_utils.get_organism_for_string(559292)
    assert result == 4932, "Failed to retrieve organism id for species id 559292."

    # --- Test case 2: Recursive lookup ---
    # We assume that taxon id 83333 is NOT directly in species data but its parent is.
    result_recursive = taxon_utils.get_organism_for_string(83333)
    assert result_recursive is None, f"Recursive lookup returned {result_recursive}, for 83333 which is not present in species data."

    # --- Test case 3: Non-existent taxon id ---
    # For an id that is very unlikely to exist, the function should return None.
    result_none = taxon_utils.get_organism_for_string(123456789)
    assert result_none is None, "Expected None for non-existent taxon id 123456789."

    # Test cases
    test_ids = [
        4932,
        9606,  # Human - should be found directly
        10090,  # Mouse - should be found directly
        511145  # E. coli K12 - should be found directly
    ]

    for test_id in test_ids:
        valid_id = taxon_utils.get_organism_for_string(test_id)
        assert valid_id == test_id, f"Expected taxon id {test_id} equal {valid_id}, but got {valid_id} for test id {test_id}."