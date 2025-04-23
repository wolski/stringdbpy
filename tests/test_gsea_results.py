import os
import tempfile
import pytest
from pathlib import Path
from loguru import logger

from string_gsea.gsea_session import GSEASession
from string_gsea.string_gsea_results import StringGSEAResults

INTEGRATION_FLAG = os.getenv("RUN_STRING_GSEA_INTEGRATION") == "1"
INTEGRATION_FLAG = True  # For testing purposes, set to True

@pytest.mark.skipif(
    not INTEGRATION_FLAG,
    reason="Integration test disabled; set RUN_STRING_GSEA_INTEGRATION=1 to enable"
)
def test_gsea_results_integration():
    # Get paths exactly as in the original code
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent
    yaml_input = project_root / "tests" / "data" / "dummy_d" / "session.yml"

    # Load session and create results object
    session = GSEASession.from_yaml(yaml_input)
    # set the base path to a new temp directory
    session.base_path = Path(tempfile.mkdtemp())
    results = StringGSEAResults(session)
    
    # Log jobs and write outputs exactly as in the original code
    logger.info(f"Jobs: {results.session.res_job_id}")
    
    # Write all outputs

    links = results.write_links()
    tsv_dir = results.write_gsea_tsv()
    graph_dir = results.write_gsea_graphs()
    
    # Add assertions to verify the test results
    for key, value in links.items():
        assert value.exists(), f"{key} file should exist"
    assert tsv_dir.exists(), "TSV directory should exist"
    assert graph_dir.exists(), "Graph directory should exist"
    
    # Check if any TSV files were created
    tsv_files = list(tsv_dir.glob("**/*.tsv"))
    assert len(tsv_files) > 0, "Should have created at least one TSV file"
    
    # Check if any graph files were created
    graph_files = list(graph_dir.glob("**/*.png"))
    assert len(graph_files) > 0, "Should have created at least one graph file"
    
    logger.info("Integration test completed successfully") 