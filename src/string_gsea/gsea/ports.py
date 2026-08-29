"""Client-owned interfaces for the GSEA application."""

from typing import Protocol

from string_gsea.gsea.model.session import CompletedJob, DownloadedJob, JobKey


class GSEAGatewayError(RuntimeError):
    """Base error raised by a GSEA gateway implementation."""


class GSEAJobFailed(GSEAGatewayError):
    """Raised when STRING-DB reports a terminal job failure."""


class GSEAGateway(Protocol):
    """Smallest STRING-DB capability exercised by the GSEA use case."""

    def submit_ranks(self, rank_data: str, species: int) -> str:
        """Submit one ranked list and return its job identifier."""
        ...

    def wait(self, job_id: str) -> CompletedJob:
        """Wait for a job to complete or raise a typed gateway error."""
        ...

    def download(self, key: JobKey, job: CompletedJob) -> DownloadedJob:
        """Download the optional artifacts exposed by a completed job."""
        ...
