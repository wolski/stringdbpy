# pyright: reportUnknownMemberType=false
"""Polymorphic rank sources and analysis policies."""

from __future__ import annotations

import io
import zipfile
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from types import MappingProxyType
from typing import Protocol

import polars as pl

from string_gsea.gsea.model.ranks import RankList, RankListCollection


class AnalysisName(StrEnum):
    """Supported differential-expression rank policies."""

    PEP_1 = "pep_1"
    PEP_1_NO_IMPUTED = "pep_1_no_imputed"
    PEP_2 = "pep_2"
    PEP_2_NO_IMPUTED = "pep_2_no_imputed"


class RankFilter(Protocol):
    """One composable transformation of differential-expression rows."""

    def apply(self, dataframe: pl.DataFrame) -> pl.DataFrame:
        """Return rows accepted by this filter."""
        ...


@dataclass(frozen=True, slots=True)
class ExcludeImputed:
    """Remove rows produced by an imputed model."""

    def apply(self, dataframe: pl.DataFrame) -> pl.DataFrame:
        return dataframe.filter(~pl.col("modelName").str.contains("(?i)imputed"))


@dataclass(frozen=True, slots=True)
class MinimumPeptides:
    """Retain rows supported by at least ``minimum`` peptides."""

    minimum: int

    def apply(self, dataframe: pl.DataFrame) -> pl.DataFrame:
        return dataframe.filter(pl.col("nrPeptides") >= self.minimum)


@dataclass(frozen=True, slots=True)
class AnalysisPolicy:
    """Named composition of rank filters."""

    name: AnalysisName
    filters: tuple[RankFilter, ...]

    def apply(self, dataframe: pl.DataFrame) -> pl.DataFrame:
        result = dataframe
        for rank_filter in self.filters:
            result = rank_filter.apply(result)
        return result


ANALYSIS_POLICIES = MappingProxyType(
    {
        AnalysisName.PEP_1: AnalysisPolicy(AnalysisName.PEP_1, ()),
        AnalysisName.PEP_1_NO_IMPUTED: AnalysisPolicy(
            AnalysisName.PEP_1_NO_IMPUTED, (ExcludeImputed(),)
        ),
        AnalysisName.PEP_2: AnalysisPolicy(AnalysisName.PEP_2, (MinimumPeptides(2),)),
        AnalysisName.PEP_2_NO_IMPUTED: AnalysisPolicy(
            AnalysisName.PEP_2_NO_IMPUTED,
            (ExcludeImputed(), MinimumPeptides(2)),
        ),
    }
)


@dataclass(frozen=True, slots=True)
class RankSourceRequest:
    """Archive plus the explicitly requested input interpretation."""

    archive: Path
    analysis: AnalysisName | None


@dataclass(frozen=True, slots=True)
class ArchiveManifest:
    """Relevant archive entries inspected once before source selection."""

    xlsx_files: tuple[str, ...]
    rank_files: tuple[str, ...]

    @classmethod
    def inspect(cls, archive: Path) -> ArchiveManifest:
        with zipfile.ZipFile(archive) as zipped:
            names = zipped.namelist()
        return cls(
            xlsx_files=tuple(name for name in names if name.endswith(".xlsx") and "DE_" in name),
            rank_files=tuple(name for name in names if name.endswith(".rnk")),
        )


class RankSource(Protocol):
    """Archive reader selected by capability rather than a mode branch."""

    def supports(self, request: RankSourceRequest, manifest: ArchiveManifest) -> bool:
        """Return whether this source can honor the request."""
        ...

    def load(self, request: RankSourceRequest, manifest: ArchiveManifest) -> RankListCollection:
        """Load ranked lists from the archive."""
        ...


class XlsxRankSource:
    """Read differential-expression XLSX data and apply an analysis policy."""

    def supports(self, request: RankSourceRequest, manifest: ArchiveManifest) -> bool:
        return request.analysis is not None and bool(manifest.xlsx_files)

    def load(self, request: RankSourceRequest, manifest: ArchiveManifest) -> RankListCollection:
        analysis = request.analysis
        if analysis is None:
            raise ValueError("XLSX rank source requires an analysis policy")
        with (
            zipfile.ZipFile(request.archive) as zipped,
            zipped.open(manifest.xlsx_files[0]) as source,
        ):
            dataframe = pl.read_excel(io.BytesIO(source.read()), sheet_name="diff_exp_analysis")
        filtered = ANALYSIS_POLICIES[analysis].apply(dataframe)
        return RankListCollection(
            analysis=analysis.value,
            rank_lists=_ranks_by_contrast(filtered),
        )


class RnkArchiveSource:
    """Read existing RNK files, including the established XLSX fallback."""

    def supports(self, request: RankSourceRequest, manifest: ArchiveManifest) -> bool:
        return bool(manifest.rank_files)

    def load(self, request: RankSourceRequest, manifest: ArchiveManifest) -> RankListCollection:
        rank_lists: list[RankList] = []
        with zipfile.ZipFile(request.archive) as zipped:
            for filename in manifest.rank_files:
                content = zipped.read(filename).decode()
                entries = {
                    fields[0]: float(fields[1])
                    for line in content.splitlines()
                    if (fields := line.split())
                }
                rank_lists.append(RankList(contrast=Path(filename).stem, entries=entries))
        return RankListCollection(analysis="from_rnk", rank_lists=rank_lists)


def select_rank_source(
    request: RankSourceRequest,
    sources: tuple[RankSource, ...] = (XlsxRankSource(), RnkArchiveSource()),
) -> RankListCollection:
    """Use the first injected source that supports the inspected archive."""
    manifest = ArchiveManifest.inspect(request.archive)
    for source in sources:
        if source.supports(request, manifest):
            return source.load(request, manifest)
    raise ValueError(f"Archive contains no supported rank input: {request.archive}")


def _ranks_by_contrast(dataframe: pl.DataFrame) -> list[RankList]:
    id_columns = [name for name in ("IDcolumn", "proteinname") if name in dataframe.columns]
    if not id_columns:
        raise ValueError("No valid ID columns found")
    identifier = id_columns[0]
    rank_lists: list[RankList] = []
    for contrast in dataframe.get_column("contrast").unique().to_list():
        rank_dataframe = dataframe.filter(pl.col("contrast") == contrast).select(
            pl.col(identifier).alias("id"), pl.col("statistic")
        )
        rank_lists.append(RankList.from_polars(rank_dataframe, contrast=str(contrast)))
    return rank_lists
