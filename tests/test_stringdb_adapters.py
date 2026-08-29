"""Tests for runtime validation at the STRING-DB HTTP boundary."""

from collections.abc import Mapping
from dataclasses import dataclass

import pytest
import requests

from string_gsea.gsea.model.session import SessionSettings
from string_gsea.gsea.ports import GSEAJobFailed
from string_gsea.ora.ports import ORAGatewayError
from string_gsea.stringdb_adapters import (
    HttpResponse,
    InvalidStringResponse,
    RequestsApiKeyProvider,
    RequestsStringDB,
    RequestsTransport,
    StringDBRequestError,
)


@dataclass(slots=True)
class FakeResponse:
    payload: object
    text: str = ""
    content: bytes = b""

    def json(self) -> object:
        return self.payload

    def raise_for_status(self) -> None:
        return None


class FakeTransport:
    def __init__(self, responses: list[HttpResponse]) -> None:
        self.responses = responses
        self.calls: list[tuple[str, str]] = []

    def _next(self) -> HttpResponse:
        return self.responses.pop(0)

    def get(self, url: str, *, params: Mapping[str, object] | None = None) -> HttpResponse:
        self.calls.append(("GET", url))
        return self._next()

    def post(self, url: str, *, data: Mapping[str, object]) -> HttpResponse:
        self.calls.append(("POST", url))
        return self._next()


def _settings() -> SessionSettings:
    return SessionSettings("secret", 0.25, 1, "tests", api_base_url="https://string.test/api")


def test_gsea_submit_poll_and_download() -> None:
    transport = FakeTransport(
        [
            FakeResponse([{"job_id": "job-1"}]),
            FakeResponse([{"status": "running"}]),
            FakeResponse(
                [
                    {
                        "status": "success",
                        "page_url": "https://page",
                        "download_url": "https://result",
                        "graph_url": "https://graph",
                    }
                ]
            ),
            FakeResponse([], text="a\tb\n"),
            FakeResponse([], content=b"PNG"),
        ]
    )
    ticks = iter((0.0, 0.0, 1.0))
    gateway = RequestsStringDB(
        _settings(),
        transport=transport,
        monotonic=lambda: next(ticks),
        sleep=lambda _seconds: None,
        poll_interval=0,
    )

    job_id = gateway.submit_ranks("P1\t1\n", 9606)
    completed = gateway.wait(job_id)
    downloaded = gateway.download(("pep_1", "A"), completed)

    assert downloaded.tsv == "a\tb\n"
    assert downloaded.graph == b"PNG"
    assert downloaded.page_url == "https://page"


@pytest.mark.parametrize(
    "payload, error",
    [
        ([], InvalidStringResponse),
        ([{"status": "error", "message": "bad input"}], GSEAJobFailed),
    ],
)
def test_submit_rejects_error_documents(payload: object, error: type[Exception]) -> None:
    gateway = RequestsStringDB(_settings(), transport=FakeTransport([FakeResponse(payload)]))
    with pytest.raises(error):
        gateway.submit_ranks("invalid", 9606)


def test_poll_failure_and_timeout() -> None:
    failed = RequestsStringDB(
        _settings(),
        transport=FakeTransport([FakeResponse([{"status": "nothing found"}])]),
        monotonic=lambda: 0,
    )
    with pytest.raises(GSEAJobFailed, match="nothing found"):
        failed.wait("job")

    timed_out = RequestsStringDB(
        _settings(),
        transport=FakeTransport([]),
        monotonic=lambda: 0,
        poll_timeout=0,
    )
    with pytest.raises(TimeoutError):
        timed_out.wait("job")


def test_taxon_and_ora_operations_validate_json() -> None:
    transport = FakeTransport(
        [
            FakeResponse([{"queryItem": "P1", "ncbiTaxonId": 9606}]),
            FakeResponse([{"queryItem": "P1", "stringId": "9606.P1"}]),
            FakeResponse([{"term": "GO:1", "genes": ["P1"]}]),
            FakeResponse(["https://network"]),
        ]
    )
    gateway = RequestsStringDB(_settings(), transport=transport)

    assert gateway.lookup_taxa(("P1",)) == {"P1": 9606}
    assert gateway.map_identifiers(("P1",), 9606) == {"P1": "9606.P1"}
    assert gateway.enrich(("9606.P1",), ("9606.P1",), 9606)[0]["term"] == "GO:1"
    assert gateway.network_link(("9606.P1",), 9606) == "https://network"


def test_ora_error_and_missing_link() -> None:
    gateway = RequestsStringDB(
        _settings(),
        transport=FakeTransport(
            [
                FakeResponse([{"error": True, "message": "bad background"}]),
                FakeResponse([]),
            ]
        ),
    )
    with pytest.raises(ORAGatewayError, match="bad background"):
        gateway.enrich(("P1",), ("P1",), 9606)
    assert gateway.network_link(("P1",), 9606) is None


@pytest.mark.parametrize(
    "payload, expected",
    [
        ([{"api_key": "key", "note": "note"}], ("key", "note")),
        ([{"api_key": "key"}], ("key", "No note provided")),
    ],
)
def test_api_key_provider_validates_and_returns_credentials(
    payload: object, expected: tuple[str, str]
) -> None:
    provider = RequestsApiKeyProvider(FakeTransport([FakeResponse(payload)]))
    assert provider.fetch() == expected


@pytest.mark.parametrize(
    "payload",
    [[], [{}], [{"api_key": "key", "note": 3}]],
)
def test_api_key_provider_rejects_invalid_payload(payload: object) -> None:
    provider = RequestsApiKeyProvider(FakeTransport([FakeResponse(payload)]))
    with pytest.raises(InvalidStringResponse):
        provider.fetch()


class FailingSession:
    def get(self, _url: str, *, params: object = None) -> requests.Response:
        raise requests.ConnectionError("offline")

    def post(self, _url: str, *, data: object) -> requests.Response:
        raise requests.ConnectionError("offline")


def test_requests_transport_translates_precise_failures(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("string_gsea.stringdb_adapters.requests.Session", FailingSession)
    transport = RequestsTransport()
    with pytest.raises(StringDBRequestError, match="offline"):
        transport.get("https://string.test")
    with pytest.raises(StringDBRequestError, match="offline"):
        transport.post("https://string.test", data={})
