"""YAML adapter preserving the established GSEA session document."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict
from pathlib import Path
from typing import cast

import yaml

from string_gsea.gsea.model.session import CompletedJob, GSEASession, JobKey, SessionSettings


class InvalidSessionDocument(ValueError):
    """Raised when a session YAML document has an invalid shape."""


def _mapping(value: object, field: str) -> Mapping[object, object]:
    if not isinstance(value, Mapping):
        raise InvalidSessionDocument(f"{field} must be a mapping")
    return cast(Mapping[object, object], value)


def _string(value: object, field: str) -> str:
    if not isinstance(value, str):
        raise InvalidSessionDocument(f"{field} must be a string")
    return value


def _number(value: object, field: str) -> int | float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise InvalidSessionDocument(f"{field} must be numeric")
    return value


def _integer(value: object, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise InvalidSessionDocument(f"{field} must be an integer")
    return value


def _decode_key(value: object) -> JobKey:
    raw = _string(value, "job key")
    parts = raw.split("~", 1)
    if len(parts) != 2:
        raise InvalidSessionDocument(f"Invalid job key: {raw}")
    return parts[0], parts[1]


def _optional_string(document: Mapping[object, object], key: str) -> str | None:
    value = document.get(key)
    if value is None:
        return None
    return _string(value, key)


def _settings(value: object) -> SessionSettings:
    document = _mapping(value, "config_dict")
    return SessionSettings(
        api_key=_string(document.get("api_key"), "api_key"),
        fdr=float(_number(document.get("fdr"), "fdr")),
        ge_enrichment_rank_direction=_integer(
            document.get("ge_enrichment_rank_direction"), "ge_enrichment_rank_direction"
        ),
        caller_identity=_string(document.get("caller_identity"), "caller_identity"),
        creation_date=_optional_string(document, "creation_date"),
        api_base_url=_string(document.get("api_base_url"), "api_base_url"),
    )


def _completed_job(value: object) -> CompletedJob:
    document = _mapping(value, "result")
    status = _string(document.get("status"), "status")
    if status != "success":
        raise InvalidSessionDocument(f"Expected completed job, got status {status}")
    return CompletedJob(
        page_url=_optional_string(document, "page_url"),
        download_url=_optional_string(document, "download_url"),
        graph_url=_optional_string(document, "graph_url"),
    )


class SessionYaml:
    """Serialize and deserialize sessions without ambiguous path-or-text arguments."""

    @staticmethod
    def dumps(session: GSEASession) -> str:
        document = {
            "current_date": session.current_date,
            "workunit_id": session.workunit_id,
            "species": session.species,
            "config_dict": asdict(session.settings),
            "base_path": str(session.base_path),
            "res_job_id": {
                f"{outer}~{inner}": job for (outer, inner), job in session.job_ids.items()
            },
            "res_data": {
                f"{outer}~{inner}": job.to_document()
                for (outer, inner), job in session.completed_jobs.items()
            },
        }
        return yaml.safe_dump(document, sort_keys=False)

    @classmethod
    def dump(cls, session: GSEASession, path: Path) -> Path:
        path.write_text(cls.dumps(session))
        return path

    @staticmethod
    def loads(content: str) -> GSEASession:
        document = _mapping(yaml.safe_load(content), "session")
        raw_jobs = _mapping(document.get("res_job_id", {}), "res_job_id")
        raw_results = _mapping(document.get("res_data", {}), "res_data")
        return GSEASession(
            current_date=_string(document.get("current_date"), "current_date"),
            workunit_id=_string(document.get("workunit_id"), "workunit_id"),
            species=_integer(document.get("species"), "species"),
            settings=_settings(document.get("config_dict")),
            base_path=Path(_string(document.get("base_path"), "base_path")),
            job_ids={_decode_key(key): _string(value, "job id") for key, value in raw_jobs.items()},
            completed_jobs={
                _decode_key(key): _completed_job(value) for key, value in raw_results.items()
            },
        )

    @classmethod
    def load(cls, path: Path) -> GSEASession:
        return cls.loads(path.read_text())
