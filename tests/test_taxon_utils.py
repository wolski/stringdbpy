"""Tests for TaxonUtils — species lookup via public API."""

import pytest

from string_gsea.taxon_utils import TaxonUtils


@pytest.fixture
def taxon_utils():
    return TaxonUtils()


def test_taxon_utils_loads_data(taxon_utils):
    """TaxonUtils should load both species and NCBI DataFrames on init."""
    assert taxon_utils.string_species_df is not None
    assert taxon_utils.string_species_df.height > 0
    assert "taxon_id" in taxon_utils.string_species_df.columns

    assert taxon_utils.ncbi_nodes_df is not None
    assert taxon_utils.ncbi_nodes_df.height > 0
    assert "taxon_id" in taxon_utils.ncbi_nodes_df.columns
    assert "parent_taxon_id" in taxon_utils.ncbi_nodes_df.columns


def test_get_organism_direct_match(taxon_utils):
    """Species directly in STRING should return themselves."""
    test_ids = [4932, 9606, 10090, 511145]
    for test_id in test_ids:
        result = taxon_utils.get_organism_for_string(test_id)
        assert result == test_id, f"Expected {test_id}, got {result}"


def test_get_organism_recursive_lookup(taxon_utils):
    """Yeast strain 559292 should resolve to parent 4932."""
    result = taxon_utils.get_organism_for_string(559292)
    assert result == 4932


def test_get_organism_not_found(taxon_utils):
    """Non-existent taxon ID should return None."""
    result = taxon_utils.get_organism_for_string(123456789)
    assert result is None


def test_get_organism_e_coli_strain(taxon_utils):
    """E. coli K12 strain 83333 — may not be in STRING, should return None or a parent."""
    result = taxon_utils.get_organism_for_string(83333)
    # 83333 is not directly in STRING species list
    assert result is None
