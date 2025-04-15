# tests/test_stringgsea.py
import pytest
from pathlib import Path

import requests
# No longer mocking requests or using os
from a373_string_gsea.stringgsea import StringGSEA

# Define the path to the test data directory relative to the test file
TEST_DATA_DIR = Path(__file__).parent / 'data'
# Corrected filename typo gesa -> gsea
TEST_JSON_PATH = TEST_DATA_DIR / 'WU_abcd_GSEA' / 'serialized_gsea_session.json'

def test_integration_deserialize_and_write_outputs(tmp_path): # Added tmp_path fixture
    """
    Integration test: Loads from an existing JSON, modifies an attribute,
    writes results and links to a temporary directory using real network calls,
    and checks if output files exist.
    """
    # Check if the real input file exists before attempting to load
    if not TEST_JSON_PATH.is_file():
        pytest.skip(f"Test skipped: Required file not found at {TEST_JSON_PATH}")

    gsea2 = StringGSEA.from_serialized_json(TEST_JSON_PATH)
    
    new_workunit_id = "xyz_integration_test" # Use a distinct ID for the test output folder
    gsea2.workunit_id = new_workunit_id

    # --- Call methods to write outputs to the temporary directory (NO MOCKING) ---
    results_output_dir = gsea2.write_gsea_results(path=tmp_path)
    links_output_files = gsea2.write_links(path=tmp_path)
    expected_base_dir = tmp_path / f"WU_{new_workunit_id}_GSEA"

    assert results_output_dir == expected_base_dir, f"write_gsea_results returned wrong base dir"
    assert expected_base_dir.is_dir(), f"Base output directory not created: {expected_base_dir}"

    # Check for each expected result file based on loaded res_data
    results_files_found = False
    for (outer_key, inner_key), data in gsea2.res_data.items():
        if data.get('status') == "success":
            # Check if download/graph URLs exist before asserting file creation
            if 'download_url' in data and 'graph_url' in data:
                results_files_found = True
                expected_sub_dir = expected_base_dir / str(outer_key)
                expected_tsv_file = expected_sub_dir / f"{inner_key}.tsv"
                expected_png_file = expected_sub_dir / f"{inner_key}.png"

                assert expected_sub_dir.is_dir(), f"Subdirectory not created: {expected_sub_dir}"
                assert expected_tsv_file.is_file(), f"Expected TSV file not found: {expected_tsv_file}"
                assert expected_tsv_file.stat().st_size > 0, f"TSV file is empty: {expected_tsv_file}" # Check if not empty
                assert expected_png_file.is_file(), f"Expected PNG file not found: {expected_png_file}"
                assert expected_png_file.stat().st_size > 0, f"PNG file is empty: {expected_png_file}" # Check if not empty
            else:
                 print(f"Skipping file check for {outer_key}/{inner_key}: Missing download_url or graph_url in loaded data.")

    if not results_files_found:
         print("Warning: No successful results with download/graph URLs found in loaded JSON data to test file writing.")
         # Depending on requirements, you might want to fail here instead:
         # pytest.fail("No successful results with download/graph URLs found in loaded JSON data.")


    # Assertions for write_links output
    links_files_found = False
    # Check if links_output_files is not empty before iterating
    if links_output_files:
        for outer_key, expected_link_file_path in links_output_files.items():
             links_files_found = True
             expected_sub_dir = expected_base_dir / str(outer_key)
             expected_file = expected_sub_dir / "links.txt"

             assert expected_link_file_path == expected_file, f"write_links returned wrong path for {outer_key}"
             assert expected_file.is_file(), f"Expected links file not found: {expected_file}"
             assert expected_file.stat().st_size > 0, f"Links file is empty: {expected_file}"
    
    # Check if get_links actually returns anything before asserting links_files_found
    if not gsea2.get_links():
         print("Warning: No links were expected based on get_links(), so no links file check performed.")
    else:
        assert links_files_found, "No links files were written (check loaded data and get_links logic)."
