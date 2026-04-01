"""Tests for StringAPIClient — HTTP interactions with STRING-DB API."""

from unittest.mock import MagicMock, patch

import pytest

from string_gsea.string_api_client import StringAPIClient


@pytest.fixture
def api_client(gsea_config):
    return StringAPIClient(gsea_config, species=9606)


@patch("string_gsea.string_api_client.requests")
def test_submit_ranks_success(mock_requests, api_client):
    mock_resp = MagicMock()
    mock_resp.json.return_value = [{"job_id": "abc123"}]
    mock_requests.post.return_value = mock_resp

    job_id = api_client.submit_ranks("GENE1\t1.5\nGENE2\t-0.3")

    assert job_id == "abc123"
    mock_requests.post.assert_called_once()


@patch("string_gsea.string_api_client.requests")
def test_submit_ranks_error(mock_requests, api_client):
    mock_resp = MagicMock()
    mock_resp.json.return_value = [{"status": "error", "message": "bad input"}]
    mock_requests.post.return_value = mock_resp

    with pytest.raises(RuntimeError, match="bad input"):
        api_client.submit_ranks("INVALID")


@patch("string_gsea.string_api_client.requests")
def test_poll_job_immediate_success(mock_requests, api_client):
    record = {"status": "success", "page_url": "http://example.com", "download_url": "http://example.com/dl"}
    mock_resp = MagicMock()
    mock_resp.json.return_value = [record]
    mock_requests.get.return_value = mock_resp

    result = api_client.poll_job("abc123", sleep_t=0, max_time=10)

    assert result["status"] == "success"
    assert result["page_url"] == "http://example.com"


@patch("string_gsea.string_api_client.requests")
def test_poll_job_failure(mock_requests, api_client):
    mock_resp = MagicMock()
    mock_resp.json.return_value = [{"status": "nothing found"}]
    mock_requests.get.return_value = mock_resp

    with pytest.raises(RuntimeError, match="nothing found"):
        api_client.poll_job("abc123", sleep_t=0, max_time=10)


@patch("string_gsea.string_api_client.requests")
def test_poll_job_timeout(mock_requests, api_client):
    mock_resp = MagicMock()
    mock_resp.json.return_value = [{"status": "running"}]
    mock_requests.get.return_value = mock_resp

    with pytest.raises(TimeoutError):
        api_client.poll_job("abc123", sleep_t=0, max_time=0)
