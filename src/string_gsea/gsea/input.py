# pyright: reportUnknownMemberType=false
"""Polymorphic rank sources and analysis policies."""

from __future__ import annotations

import io
import zipfile
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from tempfile import TemporaryDirectory
from types import MappingProxyType
from typing import Protocol

import polars as pl

from string_gsea.gsea.dea_artifact import (
    ESTIMATE_TYPE_COLUMN,
    OBSERVED_ESTIMATE,
    PEPTIDE_COUNT_COLUMN,
    DeaArtifact,
    read_dea_artifact,
)
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
class ExcludeImputedModel:
    """Remove rows whose model name marks them as imputed.

    The XLSX schema names the model per row (``Imputed_Mean_moderated``), so
    imputation is readable only from that label.
    """

    column: str = "modelName"

    def apply(self, dataframe: pl.DataFrame) -> pl.DataFrame:
        return dataframe.filter(~pl.col(self.column).str.contains("(?i)imput"))


@dataclass(frozen=True, slots=True)
class KeepObservedEstimates:
    """Keep only rows the model actually estimated from observed data.

    The AnnData schema records estimate provenance per row
    (``observed`` / ``lod_imputed`` / ``missing_fallback``), which is stricter
    and more direct than reading it off a model name.
    """

    column: str = ESTIMATE_TYPE_COLUMN
    observed: str = OBSERVED_ESTIMATE

    def apply(self, dataframe: pl.DataFrame) -> pl.DataFrame:
        if self.column not in dataframe.columns:
            raise ValueError(
                f"Cannot exclude imputed estimates: the artifact has no "
                f"{self.column!r} column, so their provenance is unknown"
            )
        return dataframe.filter(pl.col(self.column) == self.observed)


@dataclass(frozen=True, slots=True)
class MinimumPeptides:
    """Retain rows supported by at least ``minimum`` peptides."""

    minimum: int
    column: str = "nrPeptides"

    def apply(self, dataframe: pl.DataFrame) -> pl.DataFrame:
        return dataframe.filter(pl.col(self.column) >= self.minimum)


@dataclass(frozen=True, slots=True)
class RankSchema:
    """Where one input format carries the facts the policies filter on."""

    imputation: RankFilter
    peptides: str


XLSX_SCHEMA = RankSchema(imputation=ExcludeImputedModel(), peptides="nrPeptides")
ANNDATA_SCHEMA = RankSchema(imputation=KeepObservedEstimates(), peptides=PEPTIDE_COUNT_COLUMN)


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


def analysis_policy(name: AnalysisName, schema: RankSchema) -> AnalysisPolicy:
    """Compose the named policy over the columns one input format uses."""
    peptides = MinimumPeptides(2, column=schema.peptides)
    filters: dict[AnalysisName, tuple[RankFilter, ...]] = {
        AnalysisName.PEP_1: (),
        AnalysisName.PEP_1_NO_IMPUTED: (schema.imputation,),
        AnalysisName.PEP_2: (peptides,),
        AnalysisName.PEP_2_NO_IMPUTED: (schema.imputation, peptides),
    }
    return AnalysisPolicy(name, filters[name])


ANALYSIS_POLICIES = MappingProxyType(
    {name: analysis_policy(name, XLSX_SCHEMA) for name in AnalysisName}
)


@dataclass(frozen=True, slots=True)
class RankSourceRequest:
    """Archive plus the explicitly requested input interpretation."""

    archive: Path
    analysis: AnalysisName | None


@dataclass(frozen=True, slots=True)
class ArchiveManifest:
    """Relevant archive entries inspected once before source selection."""

    anndata_files: tuple[str, ...]
    xlsx_files: tuple[str, ...]
    rank_files: tuple[str, ...]

    @classmethod
    def inspect(cls, archive: Path) -> ArchiveManifest:
        with zipfile.ZipFile(archive) as zipped:
            names = zipped.namelist()
        return cls(
            anndata_files=tuple(name for name in names if name.endswith(".h5ad")),
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


class AnnDataRankSource:
    """Read the AnnData DEA artifact and rank by the column roles it records.

    Preferred over the XLSX sheet because the artifact says which column is
    the effect, the p-value and the contrast, so results from a backend with
    its own column names (SAINTexpress) rank the same way as a linear model.
    """

    def supports(self, request: RankSourceRequest, manifest: ArchiveManifest) -> bool:
        # An analysis policy is required, as for the XLSX sheet: `--which none`
        # asks for the rank files the archive already ships, not for a rederived
        # and unfiltered ranking of the same contrasts.
        return request.analysis is not None and bool(manifest.anndata_files)

    def load(self, request: RankSourceRequest, manifest: ArchiveManifest) -> RankListCollection:
        analysis = request.analysis
        if analysis is None:
            raise ValueError("AnnData rank source requires an analysis policy")
        with zipfile.ZipFile(request.archive) as zipped, TemporaryDirectory() as workdir:
            artifact = read_dea_artifact(Path(zipped.extract(manifest.anndata_files[0], workdir)))
        rows = analysis_policy(analysis, ANNDATA_SCHEMA).apply(artifact.rows)
        return RankListCollection(
            analysis=analysis.value,
            rank_lists=_ranks_from_artifact(artifact, rows),
        )


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
    sources: tuple[RankSource, ...] = (
        AnnDataRankSource(),
        XlsxRankSource(),
        RnkArchiveSource(),
    ),
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


def _ranks_from_artifact(artifact: DeaArtifact, rows: pl.DataFrame) -> list[RankList]:
    """Rank one artifact's rows exactly as prolfquapp writes its `.rnk` files.

    Identifiers that map to several subjects are averaged, which is what
    `prolfquapp:::.write_GSEA()` does before writing a rank file.
    """
    scored = (
        rows.select(
            pl.col(artifact.identifier).alias("id"),
            pl.col(artifact.roles.contrast_col).cast(pl.String).alias("contrast"),
            artifact.roles.rank_score(),
        )
        .drop_nulls()
        .filter(pl.col("score").is_finite())
    )
    averaged = scored.group_by("contrast", "id").agg(pl.col("score").mean())
    rank_lists: list[RankList] = []
    for contrast in averaged.get_column("contrast").unique().sort().to_list():
        entries = averaged.filter(pl.col("contrast") == contrast).select("id", "score")
        rank_lists.append(RankList.from_polars(entries, contrast=str(contrast)))
    return rank_lists
