import pytest
import polars as pl
from pathlib import Path
from string_gsea.gsea_utilities import get_rank_files, find_zip_files

# Define the path to the test data directory relative to the test file
TEST_DATA_DIR = Path(__file__).parent / "data"
DE_ZIP_FILENAME = "dummy_d/DEA_mouse_fasta_rnk.zip"
DE_ZIP_PATH = TEST_DATA_DIR / DE_ZIP_FILENAME


# --- Test Case 1: Test that find_zip_files raises error when no file matches pattern ---
def test_find_zip_files_raises_error_when_no_match():
    """
    Tests that find_zip_files raises FileNotFoundError when zip files exist
    in tests/data but none match the required pattern (e.g., DEA_*, not DE_*).
    """
    # Ensure the files that *don't* match the pattern exist
    non_matching_zip = TEST_DATA_DIR / "no_matching_zip" / "fasta_test.zip"
    if not non_matching_zip.exists():
        pytest.fail(
            f"Required non-matching zip file fasta_test.zip not found at {non_matching_zip}"
        )

    # Assert that the specific error is raised because no files match the DEA_* pattern
    with pytest.raises(
        FileNotFoundError, match="No zip files found matching the pattern."
    ):
        find_zip_files(TEST_DATA_DIR / "no_matching_zip")


# --- Test Case 2: Folder without zip archive ---
def test_find_zip_files_no_zip(tmp_path):
    """Tests find_zip_files raises FileNotFoundError when no zip files exist."""
    # tmp_path is an empty directory provided by pytest
    (tmp_path / "some_other_file.txt").touch()  # Add a non-zip file

    with pytest.raises(FileNotFoundError, match="No zip files found in"):
        find_zip_files(tmp_path)


# --- Test Case 3: Read rank files from the DE zip ---
def test_get_rank_files_reads_from_DE_zip():
    """Tests reading .rnk files from the DE_*.zip archive."""
    if not DE_ZIP_PATH.exists():
        pytest.fail(f"Test zip file not found at {DE_ZIP_PATH}")

    dataframes = get_rank_files(DE_ZIP_PATH)

    # Assumptions about DE_archive_with_fasta_and_rnk.zip:
    # - Contains at least one .rnk file.
    # - The .rnk file(s) are tab-separated and have 2 columns.
    # Adjust assertions based on the actual content of the test zip file.

    assert isinstance(dataframes, dict)
    assert len(dataframes) > 0, f"No .rnk files found or read from {DE_ZIP_FILENAME}."

    for key, df in dataframes.items():
        assert isinstance(key, tuple)
        assert len(key) == 2
        assert key[0] == "from_rnk"
        assert isinstance(key[1], str)  # Filename stem

        assert isinstance(df, pl.DataFrame)
        assert not df.is_empty(), (
            f"DataFrame for {key[1]}.rnk in {DE_ZIP_FILENAME} is empty."
        )
        # Assuming standard rank file format (Identifier, Rank_Value)
        assert df.shape[1] == 2, (
            f"Expected 2 columns in {key[1]}.rnk from {DE_ZIP_FILENAME}, found {df.shape[1]}"
        )
        # Add more specific content checks if needed
