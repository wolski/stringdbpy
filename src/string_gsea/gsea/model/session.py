"""Typed state produced by a completed STRING-DB GSEA run."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

type JobKey = tuple[str, str]


@dataclass(frozen=True, slots=True)
class SessionSettings:
    """Non-secret and secret settings retained by the legacy session format."""

    api_key: str
    fdr: float
    ge_enrichment_rank_direction: int
    caller_identity: str
    creation_date: str | None = None
    api_base_url: str = "https://version-12-0.string-db.org/api"


@dataclass(frozen=True, slots=True)
class CompletedJob:
    """URLs returned for a successfully completed STRING-DB job."""

    page_url: str | None = None
    download_url: str | None = None
    graph_url: str | None = None

    def to_document(self) -> dict[str, str]:
        """Return the established serialized STRING response shape."""
        document = {"status": "success"}
        for name, value in (
            ("page_url", self.page_url),
            ("download_url", self.download_url),
            ("graph_url", self.graph_url),
        ):
            if value is not None:
                document[name] = value
        return document


@dataclass(slots=True)
class GSEASession:
    """Serializable application state for one GSEA run."""

    current_date: str
    workunit_id: str
    species: int
    settings: SessionSettings
    base_path: Path
    job_ids: dict[JobKey, str] = field(default_factory=dict[JobKey, str])
    completed_jobs: dict[JobKey, CompletedJob] = field(default_factory=dict[JobKey, CompletedJob])


@dataclass(frozen=True, slots=True)
class DownloadedJob:
    """Downloaded content and links for one completed job."""

    key: JobKey
    tsv: str | None
    graph: bytes | None
    page_url: str | None


@dataclass(frozen=True, slots=True)
class GSEAArtifacts:
    """Paths and typed result returned by the GSEA application."""

    result_directory: Path
    session_path: Path
    result_json_path: Path
    archive_path: Path | None
