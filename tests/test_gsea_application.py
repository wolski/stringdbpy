"""Tests for the injected GSEA application use case."""

from datetime import datetime
from pathlib import Path
from typing import override

import pytest

from string_gsea.gsea.application import RunGSEA, RunGSEARequest
from string_gsea.gsea.model.enrichment import GSEAResult
from string_gsea.gsea.model.ranks import RankList, RankListCollection
from string_gsea.gsea.model.session import (
    CompletedJob,
    DownloadedJob,
    GSEAArtifacts,
    GSEASession,
    JobKey,
    SessionSettings,
)
from string_gsea.gsea.output import GSEAOutputWriter


class FakeGateway:
    """Deterministic GSEA gateway used as a real protocol implementation."""

    def __init__(self) -> None:
        self.submissions: list[tuple[str, int]] = []
        self.waited: list[str] = []

    def submit_ranks(self, rank_data: str, species: int) -> str:
        self.submissions.append((rank_data, species))
        return f"job-{len(self.submissions)}"

    def wait(self, job_id: str) -> CompletedJob:
        self.waited.append(job_id)
        return CompletedJob(page_url=f"https://string.test/{job_id}")

    def download(self, key: JobKey, job: CompletedJob) -> DownloadedJob:
        return DownloadedJob(key=key, tsv=None, graph=None, page_url=job.page_url)


class RecordingOutput(GSEAOutputWriter):
    """Capture application output values without serializing empty results to XLSX."""

    def __init__(self, directory: Path) -> None:
        self.directory = directory
        self.session: GSEASession | None = None
        self.result: GSEAResult | None = None

    @override
    def write(
        self,
        *,
        session: GSEASession,
        ranks: RankListCollection,
        downloads: tuple[DownloadedJob, ...],
        result: GSEAResult,
        create_zip: bool,
    ) -> GSEAArtifacts:
        self.session = session
        self.result = result
        return GSEAArtifacts(
            result_directory=self.directory,
            session_path=self.directory / "gsea_session.yml",
            result_json_path=self.directory / "result.json",
            archive_path=self.directory / "result.zip" if create_zip else None,
        )


def _settings() -> SessionSettings:
    return SessionSettings(
        api_key="secret",
        fdr=0.25,
        ge_enrichment_rank_direction=1,
        caller_identity="tests",
        api_base_url="https://string.test/api",
    )


def test_run_gsea_injects_gateway_and_returns_artifacts(tmp_path: Path) -> None:
    ranks = RankListCollection(
        "pep_1",
        [
            RankList("A", {"P1": 1.0}),
            RankList("B", {"P2": -1.0}),
        ],
    )
    gateway = FakeGateway()
    output = RecordingOutput(tmp_path)
    application = RunGSEA(
        gateway,
        output,
        now=lambda: datetime(2026, 8, 29, 12, 0, 0),
    )

    artifacts = application.execute(
        RunGSEARequest(ranks, "W1", tmp_path, 9606, _settings(), create_zip=True)
    )

    assert len(gateway.submissions) == 2
    assert gateway.waited == ["job-1", "job-2"]
    assert artifacts.archive_path == tmp_path / "result.zip"
    assert output.session is not None
    assert output.session.current_date == "2026-08-29 12:00:00"
    assert output.session.job_ids == {("pep_1", "A"): "job-1", ("pep_1", "B"): "job-2"}
    assert output.result is not None
    assert output.result.links == {
        "A_results.tsv": "https://string.test/job-1",
        "B_results.tsv": "https://string.test/job-2",
    }


class FailingGateway(FakeGateway):
    @override
    def wait(self, job_id: str) -> CompletedJob:
        raise RuntimeError(job_id)


def test_run_gsea_propagates_gateway_failure(tmp_path: Path) -> None:
    ranks = RankListCollection("pep_1", [RankList("A", {"P1": 1.0})])
    application = RunGSEA(FailingGateway(), RecordingOutput(tmp_path))

    with pytest.raises(RuntimeError, match="job-1"):
        application.execute(RunGSEARequest(ranks, "W1", tmp_path, 9606, _settings()))
