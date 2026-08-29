"""Tests for the typed TOML configuration boundary."""

from pathlib import Path

import pytest

from string_gsea import configuration
from string_gsea.configuration import GSEAConfig


def test_configuration_round_trip(tmp_path: Path) -> None:
    path = tmp_path / "config.toml"
    expected = GSEAConfig(
        api_key="secret",
        fdr=0.05,
        ge_enrichment_rank_direction=1,
        caller_identity="tests",
        creation_date="2026-08-29",
        api_base_url="https://string.test/api",
    )
    expected.write_toml(path)
    assert GSEAConfig.read_toml(path) == expected


def test_configuration_from_dict_requires_all_fields() -> None:
    with pytest.raises(ValueError, match="missing required keys"):
        GSEAConfig.from_dict({"api_key": "secret"})


@pytest.mark.parametrize(
    "field, value, message",
    [
        ("api_key", "", "api_key"),
        ("fdr", "0.25", "fdr"),
        ("ge_enrichment_rank_direction", True, "rank_direction"),
        ("caller_identity", "", "caller_identity"),
        ("creation_date", 1, "creation_date"),
        ("api_base_url", "", "api_base_url"),
    ],
)
def test_configuration_validates_boundary_types(field: str, value: object, message: str) -> None:
    data: dict[str, object] = {
        "api_key": "secret",
        "fdr": 0.25,
        "ge_enrichment_rank_direction": 1,
        "caller_identity": "tests",
    }
    data[field] = value
    with pytest.raises(ValueError, match=message):
        GSEAConfig.from_dict(data)


def test_get_configuration_reports_missing_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    missing = tmp_path / "missing.toml"
    monkeypatch.setattr(configuration, "_get_config_path", lambda: missing)
    with pytest.raises(FileNotFoundError, match="write_initial_configuration"):
        configuration.get_configuration()


def test_write_initial_configuration_writes_and_respects_declined_overwrite(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    class ApiKeys:
        def fetch(self) -> tuple[str, str]:
            return "secret", "note"

    path = tmp_path / "config.toml"
    monkeypatch.setattr(configuration, "_get_config_path", lambda: path)
    assert configuration.write_initial_configuration(ApiKeys(), "tests", 0.05) == path
    assert GSEAConfig.read_toml(path).caller_identity == "tests"

    path.write_text("unchanged")
    monkeypatch.setattr("builtins.input", lambda: "n")
    assert configuration.write_initial_configuration(ApiKeys()) == path
    assert path.read_text() == "unchanged"
