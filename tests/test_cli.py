"""Deterministic tests for root CLI composition boundaries."""

from pathlib import Path

import pytest

from string_gsea import config_cli, gsea_cli, ora_cli
from string_gsea.configuration import GSEAConfig
from string_gsea.gsea.application import RunGSEARequest
from string_gsea.gsea.input import RankSourceRequest
from string_gsea.gsea.model.ranks import RankList, RankListCollection
from string_gsea.gsea.model.session import GSEAArtifacts, SessionSettings
from string_gsea.ora.application import RunORARequest
from string_gsea.ora.model import ORAArtifacts


def _config() -> GSEAConfig:
    return GSEAConfig("secret", 0.25, 1, "tests", api_base_url="https://string.test/api")


def _fake_gateway(_settings: SessionSettings) -> object:
    return object()


class FakeSpeciesApplication:
    def __init__(self, _resolvers: object) -> None:
        pass

    def resolve(self) -> int:
        return 9606


class FakeLocalTaxonomy:
    def nearest_supported(self, taxon: int) -> int | None:
        return taxon


class FakeApiKeyProvider:
    pass


class FakeGSEAApplication:
    request_workunit: str | None = None

    def __init__(self, _gateway: object, _output: object) -> None:
        pass

    def execute(self, request: RunGSEARequest) -> GSEAArtifacts:
        FakeGSEAApplication.request_workunit = request.workunit_id
        base = Path("/tmp")
        return GSEAArtifacts(base, base / "session", base / "result", None)


def test_gsea_cli_composes_typed_application(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    archive = tmp_path / "input.zip"
    archive.touch()
    ranks = RankListCollection("from_rnk", [RankList("A", {"P1": 1.0})])

    def rank_source(_request: RankSourceRequest) -> RankListCollection:
        return ranks

    monkeypatch.setattr(gsea_cli, "get_configuration", _config)
    monkeypatch.setattr(gsea_cli, "RequestsStringDB", _fake_gateway)
    monkeypatch.setattr(gsea_cli, "select_rank_source", rank_source)
    monkeypatch.setattr(gsea_cli, "LocalTaxonomy", FakeLocalTaxonomy)
    monkeypatch.setattr(gsea_cli, "OrderedSpeciesResolver", FakeSpeciesApplication)
    monkeypatch.setattr(gsea_cli, "RunGSEA", FakeGSEAApplication)
    gsea_cli.string_gsea_run(str(archive), "W1", str(tmp_path), fdr=0.05)
    assert FakeGSEAApplication.request_workunit == "W1"


class FakeORAApplication:
    def __init__(self, _gateway: object, _output: object) -> None:
        pass

    def execute(self, request: RunORARequest) -> ORAArtifacts:
        output = request.output_directory
        return ORAArtifacts(output, {"json": output / "result.json"})


def test_ora_cli_composes_typed_application(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    significant = tmp_path / "significant.txt"
    background = tmp_path / "background.txt"
    fasta = tmp_path / "proteome.fasta"
    significant.write_text("P1\n")
    background.write_text("P1\nP2\n")
    fasta.write_text(">a OX=9606\n")
    monkeypatch.setattr(ora_cli, "get_configuration", _config)
    monkeypatch.setattr(ora_cli, "RequestsStringDB", _fake_gateway)
    monkeypatch.setattr(ora_cli, "LocalTaxonomy", FakeLocalTaxonomy)
    monkeypatch.setattr(ora_cli, "OrderedSpeciesResolver", FakeSpeciesApplication)
    monkeypatch.setattr(ora_cli, "RunORA", FakeORAApplication)
    files = ora_cli.string_ora_run(
        str(significant), str(background), str(fasta), str(tmp_path), "W1"
    )
    assert files == {"json": tmp_path / "result.json"}


def test_cli_input_guards_and_config_delegation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    with pytest.raises(FileNotFoundError, match="Zip file"):
        gsea_cli.string_gsea_run(str(tmp_path / "missing.zip"), "W1")
    with pytest.raises(FileNotFoundError, match="Significant file"):
        ora_cli.string_ora_run("missing", "background", "fasta")
    written: list[tuple[object, str, float]] = []

    def write(provider: object, identity: str, fdr: float) -> Path:
        written.append((provider, identity, fdr))
        return tmp_path / "config.toml"

    monkeypatch.setattr(config_cli, "RequestsApiKeyProvider", FakeApiKeyProvider)
    monkeypatch.setattr(config_cli, "write_initial_configuration", write)
    config_cli.write_config("tests", 0.1)
    assert len(written) == 1
    assert isinstance(written[0][0], FakeApiKeyProvider)
    assert written[0][1:] == ("tests", 0.1)
