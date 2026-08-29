"""Tests for rank-source and analysis-policy polymorphism."""

import zipfile
from pathlib import Path

import polars as pl
import pytest

from string_gsea.gsea.input import (
    ANALYSIS_POLICIES,
    AnalysisName,
    ArchiveManifest,
    RankSourceRequest,
    RnkArchiveSource,
    XlsxRankSource,
    select_rank_source,
)


def test_analysis_policies_apply_composed_filters() -> None:
    dataframe = pl.DataFrame(
        {
            "modelName": ["observed", "imputed", "observed"],
            "nrPeptides": [1, 2, 2],
        }
    )

    assert ANALYSIS_POLICIES[AnalysisName.PEP_1].apply(dataframe).height == 3
    assert ANALYSIS_POLICIES[AnalysisName.PEP_1_NO_IMPUTED].apply(dataframe).height == 2
    assert ANALYSIS_POLICIES[AnalysisName.PEP_2].apply(dataframe).height == 2
    assert ANALYSIS_POLICIES[AnalysisName.PEP_2_NO_IMPUTED].apply(dataframe).height == 1


def _rank_archive(path: Path) -> Path:
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("A.rnk", "P1\t1.5\nP2\t-0.5\n")
    return path


def test_rnk_source_and_xlsx_fallback(tmp_path: Path) -> None:
    archive = _rank_archive(tmp_path / "ranks.zip")
    ranks = select_rank_source(RankSourceRequest(archive, AnalysisName.PEP_2_NO_IMPUTED))

    assert ranks.analysis == "from_rnk"
    assert ranks["A"].entries == {"P1": 1.5, "P2": -0.5}


def test_explicit_rnk_source_reports_capability(tmp_path: Path) -> None:
    archive = _rank_archive(tmp_path / "ranks.zip")
    request = RankSourceRequest(archive, None)
    manifest = ArchiveManifest.inspect(archive)
    source = RnkArchiveSource()

    assert source.supports(request, manifest)
    assert len(source.load(request, manifest)) == 1


def test_xlsx_source_loads_real_fixture(mouse_xlsx_zip: Path) -> None:
    request = RankSourceRequest(mouse_xlsx_zip, AnalysisName.PEP_2_NO_IMPUTED)
    manifest = ArchiveManifest.inspect(mouse_xlsx_zip)
    source = XlsxRankSource()

    assert source.supports(request, manifest)
    ranks = source.load(request, manifest)
    assert ranks.analysis == AnalysisName.PEP_2_NO_IMPUTED.value
    assert len(ranks) > 0


def test_xlsx_source_requires_analysis(mouse_xlsx_zip: Path) -> None:
    request = RankSourceRequest(mouse_xlsx_zip, None)
    manifest = ArchiveManifest.inspect(mouse_xlsx_zip)

    with pytest.raises(ValueError, match="requires an analysis"):
        XlsxRankSource().load(request, manifest)


def test_no_rank_source_fails_loudly(tmp_path: Path) -> None:
    archive = tmp_path / "empty.zip"
    with zipfile.ZipFile(archive, "w"):
        pass

    with pytest.raises(ValueError, match="no supported rank input"):
        select_rank_source(RankSourceRequest(archive, None))


def test_unknown_analysis_variant_fails_at_the_boundary() -> None:
    with pytest.raises(ValueError, match="unknown"):
        AnalysisName("unknown")
