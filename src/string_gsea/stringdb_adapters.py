"""Requests-based adapters for STRING-DB external services."""

from __future__ import annotations

import time
from collections.abc import Callable, Mapping, Sequence
from typing import Protocol, cast

import requests

from string_gsea.configuration import STRING_API_BASE_DEFAULT
from string_gsea.gsea.model.session import CompletedJob, DownloadedJob, JobKey, SessionSettings
from string_gsea.gsea.ports import GSEAGatewayError, GSEAJobFailed
from string_gsea.ora.model import JsonValue, ORARecord
from string_gsea.ora.ports import ORAGatewayError

type HttpValue = str | int | float | bool


class StringDBRequestError(GSEAGatewayError, ORAGatewayError):
    """Translate a requests transport failure into client-owned gateway errors."""


class InvalidStringResponse(GSEAGatewayError, ORAGatewayError):
    """Raised when STRING-DB returns a document outside its expected schema."""


class HttpResponse(Protocol):
    """Response capability consumed by the STRING adapter."""

    @property
    def text(self) -> str: ...

    @property
    def content(self) -> bytes: ...

    def json(self) -> object: ...

    def raise_for_status(self) -> None: ...


class HttpTransport(Protocol):
    """Injectable HTTP capability used by the requests adapter."""

    def get(self, url: str, *, params: Mapping[str, HttpValue] | None = None) -> HttpResponse: ...

    def post(self, url: str, *, data: Mapping[str, HttpValue]) -> HttpResponse: ...


class RequestsTransport:
    """Production HTTP transport backed by one requests session."""

    def __init__(self) -> None:
        self._session = requests.Session()

    def get(self, url: str, *, params: Mapping[str, HttpValue] | None = None) -> requests.Response:
        try:
            return self._session.get(url, params=params)
        except requests.RequestException as error:
            raise StringDBRequestError(f"STRING request failed: {error}") from error

    def post(self, url: str, *, data: Mapping[str, HttpValue]) -> requests.Response:
        try:
            return self._session.post(url, data=data)
        except requests.RequestException as error:
            raise StringDBRequestError(f"STRING request failed: {error}") from error


def _records(payload: object) -> list[Mapping[str, object]]:
    if not isinstance(payload, list):
        raise InvalidStringResponse("STRING response must be a list")
    records: list[Mapping[str, object]] = []
    for item in cast(list[object], payload):
        if not isinstance(item, Mapping):
            raise InvalidStringResponse("STRING response entries must be mappings")
        mapping = cast(Mapping[object, object], item)
        if not all(isinstance(key, str) for key in mapping):
            raise InvalidStringResponse("STRING response entries must be string-keyed mappings")
        records.append({str(key): value for key, value in mapping.items()})
    return records


def _required_string(record: Mapping[str, object], key: str) -> str:
    value = record.get(key)
    if not isinstance(value, str) or not value:
        raise InvalidStringResponse(f"STRING response field {key} must be a non-empty string")
    return value


def _optional_string(record: Mapping[str, object], key: str) -> str | None:
    value = record.get(key)
    if value is None:
        return None
    if not isinstance(value, str):
        raise InvalidStringResponse(f"STRING response field {key} must be a string")
    return value


class RequestsApiKeyProvider:
    """Fetch configuration credentials through an injected HTTP transport."""

    def __init__(
        self,
        transport: HttpTransport | None = None,
        url: str = f"{STRING_API_BASE_DEFAULT}/json/get_api_key",
    ) -> None:
        self._transport = transport or RequestsTransport()
        self._url = url

    def fetch(self) -> tuple[str, str]:
        response = self._transport.get(self._url)
        response.raise_for_status()
        records = _records(response.json())
        if not records:
            raise InvalidStringResponse("STRING API-key response is empty")
        record = records[0]
        api_key = _required_string(record, "api_key")
        note = _optional_string(record, "note") or "No note provided"
        return api_key, note


def _json_value(value: object) -> JsonValue:
    if value is None or isinstance(value, str | int | float | bool):
        return value
    if isinstance(value, list):
        return [_json_value(item) for item in cast(list[object], value)]
    if isinstance(value, Mapping):
        mapping = cast(Mapping[object, object], value)
        if not all(isinstance(key, str) for key in mapping):
            raise InvalidStringResponse("STRING object keys must be strings")
        return {str(key): _json_value(item) for key, item in mapping.items()}
    raise InvalidStringResponse("STRING response contains a non-JSON value")


class RequestsStringDB:
    """Concrete STRING-DB GSEA and taxonomy gateway."""

    def __init__(
        self,
        settings: SessionSettings,
        *,
        sleep: Callable[[float], None] = time.sleep,
        monotonic: Callable[[], float] = time.monotonic,
        poll_interval: float = 10,
        poll_timeout: float = 3600,
        transport: HttpTransport | None = None,
    ) -> None:
        self._settings = settings
        self._sleep = sleep
        self._monotonic = monotonic
        self._poll_interval = poll_interval
        self._poll_timeout = poll_timeout
        self._transport = transport or RequestsTransport()

    def submit_ranks(self, rank_data: str, species: int) -> str:
        response = self._transport.post(
            f"{self._settings.api_base_url}/json/valuesranks_enrichment_submit",
            data={
                "species": species,
                "caller_identity": self._settings.caller_identity,
                "identifiers": rank_data,
                "api_key": self._settings.api_key,
                "ge_fdr": self._settings.fdr,
                "ge_enrichment_rank_direction": self._settings.ge_enrichment_rank_direction,
            },
        )
        response.raise_for_status()
        records = _records(response.json())
        if not records:
            raise InvalidStringResponse("STRING submit response is empty")
        record = records[0]
        if record.get("status") == "error":
            raise GSEAJobFailed(_optional_string(record, "message") or "STRING submit failed")
        return _required_string(record, "job_id")

    def wait(self, job_id: str) -> CompletedJob:
        started = self._monotonic()
        while self._monotonic() - started < self._poll_timeout:
            response = self._transport.get(
                f"{self._settings.api_base_url}/json/valuesranks_enrichment_status",
                params={"api_key": self._settings.api_key, "job_id": job_id},
            )
            response.raise_for_status()
            records = _records(response.json())
            if not records:
                raise InvalidStringResponse("STRING poll response is empty")
            record = records[0]
            status = _required_string(record, "status")
            if status == "success":
                return CompletedJob(
                    page_url=_optional_string(record, "page_url"),
                    download_url=_optional_string(record, "download_url"),
                    graph_url=_optional_string(record, "graph_url"),
                )
            if status in {"nothing found", "unknown organism", "error"}:
                raise GSEAJobFailed(f"STRING job {job_id} failed: {status}")
            self._sleep(self._poll_interval)
        raise TimeoutError(f"STRING job {job_id} did not finish in {self._poll_timeout:g}s")

    def download(self, key: JobKey, job: CompletedJob) -> DownloadedJob:
        tsv = self._download_text(job.download_url)
        graph = self._download_bytes(job.graph_url)
        return DownloadedJob(key=key, tsv=tsv, graph=graph, page_url=job.page_url)

    def lookup_taxa(self, identifiers: Sequence[str]) -> dict[str, int | None]:
        response = self._transport.get(
            f"{self._settings.api_base_url}/json/get_string_ids",
            params={"identifiers": "\r".join(identifiers)},
        )
        response.raise_for_status()
        result: dict[str, int | None] = {}
        for record in _records(response.json()):
            identifier = record.get("queryItem") or record.get("inputId")
            taxon = record.get("ncbiTaxonId")
            if isinstance(identifier, str) and (isinstance(taxon, int) or taxon is None):
                result[identifier] = taxon
        return result

    def map_identifiers(self, identifiers: Sequence[str], species: int) -> dict[str, str]:
        result: dict[str, str] = {}
        for offset in range(0, len(identifiers), 500):
            batch = identifiers[offset : offset + 500]
            response = self._transport.get(
                f"{self._settings.api_base_url}/json/get_string_ids",
                params={
                    "identifiers": "\r".join(batch),
                    "species": species,
                    "caller_identity": self._settings.caller_identity,
                },
            )
            response.raise_for_status()
            for record in _records(response.json()):
                query = record.get("queryItem")
                string_id = record.get("stringId")
                if isinstance(query, str) and isinstance(string_id, str):
                    result[query] = string_id
        return result

    def enrich(
        self,
        identifiers: Sequence[str],
        background: Sequence[str],
        species: int,
    ) -> tuple[ORARecord, ...]:
        response = self._transport.post(
            f"{self._settings.api_base_url}/json/enrichment",
            data={
                "identifiers": "\r".join(identifiers),
                "background_string_identifiers": "\r".join(background),
                "species": species,
                "caller_identity": self._settings.caller_identity,
            },
        )
        response.raise_for_status()
        records = _records(response.json())
        if records and records[0].get("error"):
            raise ORAGatewayError(
                _optional_string(records[0], "message") or "STRING enrichment failed"
            )
        return tuple(
            {key: _json_value(value) for key, value in record.items()} for record in records
        )

    def network_link(self, identifiers: Sequence[str], species: int) -> str | None:
        response = self._transport.post(
            f"{self._settings.api_base_url}/json/get_link",
            data={
                "identifiers": "\r".join(identifiers),
                "species": species,
                "caller_identity": self._settings.caller_identity,
            },
        )
        response.raise_for_status()
        payload: object = response.json()
        if not isinstance(payload, list) or not payload:
            return None
        link: object = cast(list[object], payload)[0]
        if not isinstance(link, str):
            raise InvalidStringResponse("STRING network link must be a string")
        return link

    def _download_text(self, url: str | None) -> str | None:
        if url is None:
            return None
        response = self._transport.get(url)
        response.raise_for_status()
        return response.text

    def _download_bytes(self, url: str | None) -> bytes | None:
        if url is None:
            return None
        response = self._transport.get(url)
        response.raise_for_status()
        return response.content
