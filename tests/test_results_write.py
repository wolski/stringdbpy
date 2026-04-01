"""Tests for StringGSEAResults — download and write methods with mocked HTTP."""

from unittest.mock import MagicMock, patch

import pytest

from string_gsea.gsea_session import GSEASession
from string_gsea.string_gsea_results import StringGSEAResults


@pytest.fixture
def fake_session(tmp_path, gsea_config):
    """GSEASession with fake res_data containing mock URLs."""
    return GSEASession(
        current_date="2026-04-01",
        workunit_id="TEST01",
        species=9606,
        config_dict=gsea_config,
        base_path=tmp_path,
        res_job_id={("pep_1", "contrast_A"): "job_1"},
        res_data={
            ("pep_1", "contrast_A"): {
                "status": "success",
                "page_url": "http://string-db.org/page/123",
                "download_url": "http://string-db.org/dl/results.tsv",
                "graph_url": "http://string-db.org/dl/graph.png",
            }
        },
    )


@pytest.fixture
def results(fake_session):
    return StringGSEAResults(fake_session)


def test_get_links_extracts_urls(results):
    links = results.get_links()
    assert "pep_1" in links
    assert "contrast_A" in links["pep_1"]
    assert links["pep_1"]["contrast_A"] == "http://string-db.org/page/123"


def test_write_links_creates_files(results, tmp_path):
    results.write_links(tmp_path)

    links_file = tmp_path / "pep_1" / "links.txt"
    assert links_file.exists()
    content = links_file.read_text()
    assert "contrast_A" in content
    assert "http://string-db.org/page/123" in content


@patch("string_gsea.string_gsea_results.requests")
def test_download_populates_tsv_and_graph_content(mock_requests, results):
    mock_tsv_resp = MagicMock()
    mock_tsv_resp.text = "col1\tcol2\nval1\tval2\n"
    mock_graph_resp = MagicMock()
    mock_graph_resp.content = b"\x89PNG fake image data"
    mock_requests.get.side_effect = [mock_tsv_resp, mock_graph_resp]

    results.download()

    assert ("pep_1", "contrast_A") in results.tsv_content
    assert results.tsv_content[("pep_1", "contrast_A")] == "col1\tcol2\nval1\tval2\n"
    assert ("pep_1", "contrast_A") in results.graph_content
    assert results.graph_content[("pep_1", "contrast_A")] == b"\x89PNG fake image data"


def test_write_tsv_writes_from_memory(results, tmp_path):
    results.tsv_content = {("pep_1", "contrast_A"): "col1\tcol2\nval1\tval2\n"}

    results.write_tsv(tmp_path)

    tsv_files = list(tmp_path.glob("**/*_results.tsv"))
    assert len(tsv_files) == 1
    assert tsv_files[0].read_text() == "col1\tcol2\nval1\tval2\n"


def test_write_graphs_writes_from_memory(results, tmp_path):
    results.graph_content = {("pep_1", "contrast_A"): b"\x89PNG fake image data"}

    results.write_graphs(tmp_path)

    png_files = list(tmp_path.glob("**/*_results.png"))
    assert len(png_files) == 1
    assert png_files[0].read_bytes() == b"\x89PNG fake image data"


def test_zip_folder(tmp_path):
    # Create a folder with some files
    folder = tmp_path / "test_folder"
    folder.mkdir()
    (folder / "file1.txt").write_text("hello")
    (folder / "file2.txt").write_text("world")

    archive_path = StringGSEAResults.zip_folder(folder)

    assert archive_path.exists()
    assert archive_path.suffix == ".zip"
