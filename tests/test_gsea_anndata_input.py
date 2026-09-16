"""Tests for the AnnData DEA-artifact rank source, across every prolfqua model.

Each fixture archive was produced by one prolfqua modelling facade (see
`scripts/generate_h5ad_fixtures.R`) and carries both the `AnnData.h5ad` this
source reads and the `.rnk` files prolfquapp derived from the same contrasts.
Comparing the two is what keeps one definition of "the rank" across R and
Python.
"""

import math
import re
import zipfile
from pathlib import Path

import pytest

from string_gsea.gsea.dea_artifact import read_dea_artifact
from string_gsea.gsea.input import (
    AnalysisName,
    AnnDataRankSource,
    ArchiveManifest,
    RankSourceRequest,
    select_rank_source,
)
from string_gsea.gsea.model.ranks import RankListCollection
from tests.conftest import anndata_fixture_archives

ARCHIVES = anndata_fixture_archives()
MODELS = [archive.stem for archive in ARCHIVES]

# prolfquapp names a rank file `GSEA_<contrast>_WU<workunit>.rnk`, or
# `Bait_<contrast>.rnk` for a directional backend.
RNK_NAME = re.compile(r"^(?:GSEA_|Bait_)(?P<contrast>.+?)(?:_WU.*)?$")


def _load(archive: Path, analysis: AnalysisName | None = None) -> RankListCollection:
    request = RankSourceRequest(archive, analysis)
    manifest = ArchiveManifest.inspect(archive)
    return AnnDataRankSource().load(request, manifest)


def _rnk_ranks(archive: Path) -> dict[str, dict[str, float]]:
    """The rank files prolfquapp wrote into the same archive."""
    ranks: dict[str, dict[str, float]] = {}
    with zipfile.ZipFile(archive) as zipped:
        for name in zipped.namelist():
            if not name.endswith(".rnk"):
                continue
            match = RNK_NAME.match(Path(name).stem)
            assert match is not None, name
            ranks[match.group("contrast")] = {
                fields[0]: float(fields[1])
                for line in zipped.read(name).decode().splitlines()
                if (fields := line.split())
            }
    return ranks


def test_every_facade_has_a_fixture() -> None:
    assert len(ARCHIVES) >= 20, "regenerate with scripts/generate_h5ad_fixtures.R"
    assert "saint" in MODELS
    assert "lm" in MODELS


@pytest.mark.parametrize("archive", ARCHIVES, ids=MODELS)
def test_anndata_source_ranks_every_model(archive: Path) -> None:
    ranks = _load(archive)

    assert ranks.analysis == "from_anndata"
    assert len(ranks) > 0
    for rank_list in ranks:
        assert rank_list.n_genes > 0
        assert all(isinstance(score, float) for score in rank_list.entries.values())


@pytest.mark.parametrize("archive", ARCHIVES, ids=MODELS)
def test_anndata_ranks_match_the_rank_files_prolfquapp_wrote(archive: Path) -> None:
    ranks = _load(archive)
    expected = _rnk_ranks(archive)

    assert set(ranks.contrasts) == set(expected)
    for contrast, entries in expected.items():
        produced = ranks[contrast].entries
        assert set(produced) == set(entries)
        for identifier, score in entries.items():
            assert math.isclose(produced[identifier], score, rel_tol=1e-9, abs_tol=1e-12)


@pytest.mark.parametrize("archive", ARCHIVES, ids=MODELS)
def test_anndata_source_wins_over_the_other_inputs(archive: Path) -> None:
    # The fixtures also carry .rnk files, so this pins the selection order.
    ranks = select_rank_source(RankSourceRequest(archive, None))

    assert ranks.analysis == "from_anndata"


def test_saint_ranks_on_the_effect_size_it_reports(tmp_path: Path) -> None:
    archive = next(a for a in ARCHIVES if a.stem == "saint")
    with zipfile.ZipFile(archive) as zipped:
        extracted = Path(zipped.extract("AnnData.h5ad", path=str(tmp_path)))
    artifact = read_dea_artifact(extracted)

    roles = artifact.roles
    assert roles.contrast_col == "Bait"
    assert roles.effect_col == "log2_EFCs"
    assert roles.score_col == "SaintScore"
    assert roles.pvalue_col is None
    assert not roles.tests_differences
    assert roles.directional
    # No p-value means the effect size is the rank, not the bounded score.
    ranks = _load(archive)
    effects = artifact.rows.select("log2_EFCs").to_series().to_list()
    assert math.isclose(min(ranks["A"].entries.values()), min(effects), rel_tol=1e-9)


def test_linear_model_ranks_are_signed_log_p_values() -> None:
    archive = next(a for a in ARCHIVES if a.stem == "lm")
    ranks = _load(archive)

    assert ranks.analysis == "from_anndata"
    scores = [score for rank_list in ranks for score in rank_list.entries.values()]
    # Signed -log10(p) straddles zero and leaves the |p| = 1 features at zero.
    assert min(scores) < 0 < max(scores)


@pytest.mark.parametrize(
    "model",
    ["lm_impute", "limma_impute", "rfit_impute", "limma_voom_impute", "lm_missing"],
)
def test_no_imputed_policy_drops_filled_in_estimates(model: str) -> None:
    archive = next(a for a in ARCHIVES if a.stem == model)

    everything = _load(archive, AnalysisName.PEP_1)
    observed_only = _load(archive, AnalysisName.PEP_1_NO_IMPUTED)

    assert observed_only.analysis == AnalysisName.PEP_1_NO_IMPUTED.value
    total = sum(rank_list.n_genes for rank_list in everything)
    kept = sum(rank_list.n_genes for rank_list in observed_only)
    assert 0 < kept < total


@pytest.mark.parametrize("archive", ARCHIVES, ids=MODELS)
def test_minimum_peptides_policy_reads_the_annotation(archive: Path) -> None:
    every_peptide = _load(archive, AnalysisName.PEP_1)
    two_peptides = _load(archive, AnalysisName.PEP_2)

    total = sum(rank_list.n_genes for rank_list in every_peptide)
    kept = sum(rank_list.n_genes for rank_list in two_peptides)
    assert 0 < kept <= total


def test_reading_a_non_prolfquapp_anndata_fails_loudly(tmp_path: Path) -> None:
    anndata = pytest.importorskip("anndata")
    import numpy as np

    path = tmp_path / "bare.h5ad"
    anndata.AnnData(X=np.zeros((2, 2), dtype=float)).write_h5ad(path)

    with pytest.raises(ValueError, match="not written by prolfquapp"):
        read_dea_artifact(path)
