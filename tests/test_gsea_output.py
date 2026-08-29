"""Compatibility tests for concrete GSEA filesystem output."""

from pathlib import Path

import pytest

from string_gsea.gsea.model.enrichment import GSEAResult
from string_gsea.gsea.model.ranks import RankList, RankListCollection
from string_gsea.gsea.model.session import DownloadedJob, GSEASession, SessionSettings
from string_gsea.gsea.output import GSEAOutputWriter

COMPATIBILITY = Path(__file__).parent / "data" / "compatibility"


def test_output_writer_preserves_directory_and_download_layout(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    ranks = RankListCollection("pep_1", [RankList("A", {"P1": 1.0})])
    result = GSEAResult(data={}, rank_lists={"A": ranks.first()})
    session = GSEASession(
        current_date="2026-08-29 12:00:00",
        workunit_id="W1",
        species=9606,
        settings=SessionSettings("secret", 0.25, 1, "tests"),
        base_path=tmp_path,
    )
    xlsx_calls: list[tuple[str, Path]] = []

    def write_xlsx(_result: GSEAResult, workunit_id: str, output: Path) -> None:
        xlsx_calls.append((workunit_id, output))

    monkeypatch.setattr("string_gsea.gsea.output.write_gsea_xlsx", write_xlsx)
    artifacts = GSEAOutputWriter().write(
        session=session,
        ranks=ranks,
        downloads=(
            DownloadedJob(
                ("pep_1", "A"),
                "header\nvalue\n",
                b"png",
                "https://string.test/A",
            ),
            DownloadedJob(("pep_1", "B"), None, None, None),
        ),
        result=result,
        create_zip=True,
    )

    output = tmp_path / "WU_W1_GSEA"
    assert artifacts.result_directory == output
    assert artifacts.archive_path == tmp_path / "WU_W1_GSEA.zip"
    assert artifacts.archive_path is not None
    assert artifacts.archive_path.exists()
    assert artifacts.session_path.name == "gsea_session.yml"
    assert artifacts.result_json_path.name == "WUW1_gsea_result.json"
    assert artifacts.result_json_path.read_text() == (
        COMPATIBILITY / "gsea_result.json"
    ).read_text().rstrip("\n")
    assert (output / "pep_1" / "A.rnk").read_text() == "P1\t1.0\n"
    assert (output / "pep_1" / "A_results.tsv").read_text() == "header\nvalue\n"
    assert (output / "pep_1" / "A_results.png").read_bytes() == b"png"
    assert (output / "pep_1" / "links.txt").read_text() == "A: https://string.test/A\n"
    assert xlsx_calls == [("W1", output)]
    output_paths = {
        "result_directory": artifacts.result_directory.relative_to(tmp_path).as_posix(),
        "session": artifacts.session_path.relative_to(tmp_path).as_posix(),
        "json": artifacts.result_json_path.relative_to(tmp_path).as_posix(),
        "archive": artifacts.archive_path.relative_to(tmp_path).as_posix(),
    }
    assert output_paths == {
        "result_directory": "WU_W1_GSEA",
        "session": "WU_W1_GSEA/gsea_session.yml",
        "json": "WU_W1_GSEA/WUW1_gsea_result.json",
        "archive": "WU_W1_GSEA.zip",
    }
