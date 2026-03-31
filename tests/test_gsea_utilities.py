import polars as pl
import pytest

from string_gsea.gsea_utilities import get_rank_files


# --- Test Case: Read rank files from the DE zip ---
def test_get_rank_files_reads_from_DE_zip(dummy_session_rnk_zip):
    """Tests reading .rnk files from the DE_*.zip archive."""
    if not dummy_session_rnk_zip.exists():
        pytest.fail(f"Test zip file not found at {dummy_session_rnk_zip}")

    dataframes = get_rank_files(dummy_session_rnk_zip)

    # Assumptions about DE_archive_with_fasta_and_rnk.zip:
    # - Contains at least one .rnk file.
    # - The .rnk file(s) are tab-separated and have 2 columns.
    # Adjust assertions based on the actual content of the test zip file.

    assert isinstance(dataframes, dict)
    assert len(dataframes) > 0, f"No .rnk files found or read from {dummy_session_rnk_zip}."

    for key, df in dataframes.items():
        assert isinstance(key, tuple)
        assert len(key) == 2
        assert key[0] == "from_rnk"
        assert isinstance(key[1], str)  # Filename stem

        assert isinstance(df, pl.DataFrame)
        assert not df.is_empty(), f"DataFrame for {key[1]}.rnk is empty."
        # Assuming standard rank file format (Identifier, Rank_Value)
        assert df.shape[1] == 2, f"Expected 2 columns in {key[1]}.rnk, found {df.shape[1]}"
        # Add more specific content checks if needed
