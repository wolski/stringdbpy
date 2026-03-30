
import polars as pl
import pytest

from string_gsea import network


@pytest.fixture
def long_results_df(network_test_xlsx) -> pl.DataFrame:
    # Load a small sample of the Excel file for testing
    df = pl.read_excel(network_test_xlsx)
    cont = "Bait_NCP_pUbT12T14_results.tsv"
    cat = "Monarch"
    df = df.filter((pl.col("contrast") == cont) & (pl.col("category") == cat))
    return df


@pytest.fixture
def processed_df(long_results_df: pl.DataFrame) -> pl.DataFrame:
    """DataFrame after explode_protein_columns."""
    return network.explode_protein_columns(long_results_df)


@pytest.fixture
def summarized_df(processed_df: pl.DataFrame) -> pl.DataFrame:
    """DataFrame after summarize_terms."""
    return network.summarize_terms(processed_df)


def test_select_top_terms(long_results_df: pl.DataFrame):
    """Test that select_top_terms limits rows per (contrast, category)."""
    result = network.select_top_terms(long_results_df, max_terms=3)
    assert isinstance(result, pl.DataFrame)
    assert result.height <= 3
    # Verify ordering by FDR (first row should have lowest FDR)
    fdrs = result["falseDiscoveryRate"].to_list()
    assert fdrs == sorted(fdrs)


def test_select_top_terms_passthrough(long_results_df: pl.DataFrame):
    """Test that max_terms larger than data returns all rows."""
    result = network.select_top_terms(long_results_df, max_terms=10000)
    assert result.height == long_results_df.height


def test_explode_protein_columns(
    processed_df: pl.DataFrame, long_results_df: pl.DataFrame
):
    assert isinstance(processed_df, pl.DataFrame)
    # Check that columns are exploded to lists
    assert "proteinIDs" in processed_df.columns
    assert processed_df.height > long_results_df.height


def test_summarize_terms(summarized_df: pl.DataFrame):
    assert isinstance(summarized_df, pl.DataFrame)
    # Should have filtered for FDR and genesMapped
    assert (summarized_df["falseDiscoveryRate"] < 0.05).all()
    assert (summarized_df["genesMapped"] > 10).all()
