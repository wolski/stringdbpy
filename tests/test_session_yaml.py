"""Compatibility and validation tests for the session YAML adapter."""

from pathlib import Path

import pytest

from string_gsea.gsea.model.session import CompletedJob, GSEASession, SessionSettings
from string_gsea.gsea.session_yaml import InvalidSessionDocument, SessionYaml

COMPATIBILITY = Path(__file__).parent / "data" / "compatibility" / "session.yml"


def _session() -> GSEASession:
    return GSEASession(
        current_date="2026-08-29 12:00:00",
        workunit_id="W1",
        species=9606,
        settings=SessionSettings(
            api_key="secret",
            fdr=0.25,
            ge_enrichment_rank_direction=1,
            caller_identity="tests",
            api_base_url="https://string.test/api",
        ),
        base_path=Path("/tmp/string-gsea"),
        job_ids={("pep_1", "A"): "job-1"},
        completed_jobs={
            ("pep_1", "A"): CompletedJob(
                page_url="https://string.test/page",
                download_url="https://string.test/result.tsv",
                graph_url="https://string.test/graph.png",
            )
        },
    )


def test_session_yaml_is_byte_compatible() -> None:
    assert SessionYaml.dumps(_session()) == COMPATIBILITY.read_text()


def test_session_yaml_has_explicit_text_and_path_operations(tmp_path: Path) -> None:
    path = SessionYaml.dump(_session(), tmp_path / "session.yml")

    assert SessionYaml.loads(path.read_text()) == _session()
    assert SessionYaml.load(path) == _session()


@pytest.mark.parametrize(
    "content, message",
    [
        ("[]\n", "session must be a mapping"),
        ("current_date: now\n", "workunit_id must be a string"),
        (COMPATIBILITY.read_text().replace("species: 9606", "species: true"), "species"),
        (COMPATIBILITY.read_text().replace("status: success", "status: failed"), "completed job"),
    ],
)
def test_session_yaml_rejects_invalid_documents(content: str, message: str) -> None:
    with pytest.raises(InvalidSessionDocument, match=message):
        SessionYaml.loads(content)
