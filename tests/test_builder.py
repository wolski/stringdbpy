"""Tests for StringGSEABuilder — orchestration with mocked API client."""

from unittest.mock import MagicMock

import pytest

from string_gsea.gsea_session import GSEASession
from string_gsea.models.gsea_models import RankList, RankListCollection
from string_gsea.string_gsea_builder import StringGSEABuilder


@pytest.fixture
def rank_lists():
    """Two small RankLists in a RankListCollection."""
    rl1 = RankList(contrast="contrast_A", entries={"GENE1": 1.5, "GENE2": -0.3})
    rl2 = RankList(contrast="contrast_B", entries={"GENE3": 2.1, "GENE4": 0.7})
    return RankListCollection(analysis="pep_1", rank_lists=[rl1, rl2])


@pytest.fixture
def builder(rank_lists, gsea_config, tmp_path):
    return StringGSEABuilder(
        rank_lists=rank_lists,
        config=gsea_config,
        workunit_id="TEST01",
        species=9606,
        base_path=tmp_path,
    )


def test_submit_calls_api_for_each_rank_list(builder):
    builder.api_client = MagicMock()
    builder.api_client.submit_ranks.side_effect = ["job_1", "job_2"]

    builder.submit()

    assert builder.api_client.submit_ranks.call_count == 2
    assert len(builder.session.res_job_id) == 2
    assert set(builder.session.res_job_id.values()) == {"job_1", "job_2"}


def test_submit_sends_rnk_string(builder):
    builder.api_client = MagicMock()
    builder.api_client.submit_ranks.return_value = "job_1"

    builder.submit()

    # Check that submit_ranks received a tab-separated string
    call_args = builder.api_client.submit_ranks.call_args_list[0]
    rank_str = call_args[0][0]
    assert "\t" in rank_str
    assert rank_str.endswith("\n")


def test_poll_populates_res_data(builder):
    builder.session.res_job_id = {("pep_1", "contrast_A"): "job_1"}
    mock_result = {"status": "success", "page_url": "http://example.com"}
    builder.api_client = MagicMock()
    builder.api_client.poll_job.return_value = mock_result

    builder.poll()

    assert builder.session.res_data[("pep_1", "contrast_A")] == mock_result


def test_poll_without_submit_raises(builder):
    with pytest.raises(RuntimeError, match="No jobs to poll"):
        builder.poll()


def test_submit_poll_chain_returns_builder(builder):
    builder.api_client = MagicMock()
    builder.api_client.submit_ranks.side_effect = ["job_1", "job_2"]
    builder.api_client.poll_job.return_value = {"status": "success", "page_url": "http://example.com"}

    result = builder.submit().poll()
    assert isinstance(result, StringGSEABuilder)
    assert len(builder.session.res_data) == 2


def test_save_session(builder, tmp_path):
    builder.session.res_job_id = {("pep_1", "contrast_A"): "job_1"}
    builder.session.res_data = {("pep_1", "contrast_A"): {"status": "success"}}

    yaml_path = builder.save_session(tmp_path / "gsea_session.yml")

    assert yaml_path.exists()
    loaded = GSEASession.from_yaml(yaml_path)
    assert loaded.workunit_id == "TEST01"
    assert ("pep_1", "contrast_A") in loaded.res_job_id
