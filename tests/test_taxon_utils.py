import pytest
import polars as pl
from pathlib import Path

# Import functions and relevant file paths from your module.
# Adjust the import if your module name is different.
from src.a373_string_gsea.taxon_utils import (
    read_species_string_data,
    read_ncbi_nodes_data,
    get_organism_for_string,
    SPECIES_ZIP_PATH,
    NCBI_ZIP_PATH
)

# Check file existence for skipping tests if necessary.
species_file_exists = SPECIES_ZIP_PATH.exists()
ncbi_file_exists = NCBI_ZIP_PATH.exists()


@pytest.mark.skipif(not species_file_exists, reason="Species zip file not found")
def test_read_species_string_data():
    species_df = read_species_string_data()
    assert species_df is not None, "The species DataFrame should not be None."

    # Check that the '#taxon_id' column was renamed to "taxon_id"
    assert "taxon_id" in species_df.columns, "Expected column 'taxon_id' was not found."
    assert "#taxon_id" not in species_df.columns, "Column '#taxon_id' was not renamed."

    # Optionally, ensure that some rows are loaded.
    assert species_df.height > 0, "No rows found in species data."


@pytest.mark.skipif(not ncbi_file_exists, reason="NCBI zip file not found")
def test_read_ncbi_nodes_data():
    ncbi_df = read_ncbi_nodes_data()
    assert ncbi_df is not None, "The NCBI DataFrame should not be None."

    # Check that both required columns are present.
    for col in ("taxon_id", "parent_taxon_id"):
        assert col in ncbi_df.columns, f"Expected column '{col}' not found in NCBI data."

    # Optionally, check that the file has some rows.
    assert ncbi_df.height > 0, "No rows found in NCBI data."


@pytest.mark.skipif(not (species_file_exists and ncbi_file_exists), reason="Necessary data files not found")
def test_get_species_for_string():
    species_df = read_species_string_data()
    ncbi_df = read_ncbi_nodes_data()

    # Convert species taxon ids to a Python list.
    species_ids = species_df["taxon_id"].to_list()

    # --- Test case 1: Direct match ---
    # Here, we assume that taxon id 559292 is present in species data.
    result = get_organism_for_string(559292, ncbi_df, species_df)
    assert result == 4932, "Failed to retrieve organism id for species id 559292."

    # --- Test case 2: Recursive lookup ---
    # We assume that taxon id 83333 is NOT directly in species data but its parent is.
    result_recursive = get_organism_for_string(83333, ncbi_df, species_df)
    assert result_recursive is None, f"Recursive lookup returned {result_recursive}, for 83333 which is not present in species data."

    # --- Test case 3: Non-existent taxon id ---
    # For an id that is very unlikely to exist, the function should return None.
    result_none = get_organism_for_string(123456789, ncbi_df, species_df)
    assert result_none is None, "Expected None for non-existent taxon id 123456789."

