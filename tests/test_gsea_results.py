import tempfile
from pathlib import Path

import pytest
from loguru import logger

from string_gsea.gsea_session import GSEASession
from string_gsea.string_gsea_results import StringGSEAResults


@pytest.mark.integration
def test_gsea_results_integration(dummy_session_yml):
    # Load session and create results object
    session = GSEASession.from_yaml(dummy_session_yml)
    # set the base path to a new temp directory
    session.base_path = Path(tempfile.mkdtemp())
    results = StringGSEAResults(session)

    # Log jobs and write outputs exactly as in the original code
    logger.info(f"Jobs: {results.session.res_job_id}")

    # Download into memory
    results.download()

    # Write all outputs
    out_dir = session.base_path / f"WU_{session.workunit_id}_GSEA"
    out_dir.mkdir(parents=True, exist_ok=True)

    results.write_links(out_dir)
    results.write_tsv(out_dir)
    results.write_graphs(out_dir)

    # Check if any TSV files were created
    tsv_files = list(out_dir.glob("**/*.tsv"))
    assert len(tsv_files) > 0, "Should have created at least one TSV file"

    # Check if any graph files were created
    graph_files = list(out_dir.glob("**/*.png"))
    assert len(graph_files) > 0, "Should have created at least one graph file"

    logger.info("Integration test completed successfully")
