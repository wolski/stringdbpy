# pyright: reportUnknownMemberType=false, reportUnknownVariableType=false
# pyright: reportUnknownArgumentType=false, reportMissingTypeStubs=false
"""Reading a prolfquapp differential-expression result artifact (AnnData).

prolfquapp writes `AnnData.h5ad` beside its report. The file is self-describing:
`uns["prolfquapp"]["contrast_configuration"]` records which column plays which
role for the modelling backend that ran, so nothing here needs to know whether
the numbers came from a linear model, limma or SAINTexpress.

Layout written by `prolfquapp::as_AnnData.SummarizedExperiment()`:

- `var` carries the per-feature annotation (identifiers, peptide counts)
- each contrast is one `varm` entry, listed in order by `varm_key_order`
- `uns["prolfquapp"]` carries the column roles, provenance and schema version

The `varm` entry has two layouts. prolfquapp 2.10.x writes a numeric matrix,
its columns named in `uns` by `varm_columns` and its non-numeric columns in
`varm_annotations`. From 2.11.0 it writes a data frame carrying its own columns,
the feature keys among them, and `varm_columns` is gone.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol, cast

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
from numpy.typing import NDArray

# The annotation column holding the identifier enrichment tools are given.
# Which column that is depends on the reader that produced the analysis
# (`IDcolumn`, `cleanID`, ...), so prolfquapp records the name in the artifact.
IDENTIFIER_KEY = "identifier_key"

# Row-level estimate provenance: rows that are not observed were filled in by
# the model (imputation, or a fallback for an unestimable contrast).
ESTIMATE_TYPE_COLUMN = "estimate_type"
OBSERVED_ESTIMATE = "observed"

# Per-feature peptide count, as prolfquapp names it in the annotation frame.
# prolfquapp 2.10.4 settled on `nrPeptides` everywhere; analyses written before
# it spelled the same count `nr_peptides` on the simulated and MSstats reader
# paths. Both names are read, in this order, because results already published
# cannot be rewritten.
PEPTIDE_COUNT_COLUMNS = ("nrPeptides", "nr_peptides")

_PROLFQUAPP_UNS = "prolfquapp"
_DEA_ARTIFACT_TYPE = "dea_results"

# AnnData has no missing value for strings. R's HDF5 writer stores one as the
# literal "NA" plus an `rhdf5-NA.OK` attribute, which `anndata` does not act
# on, so a missing identifier arrives here as that text.
_R_MISSING_STRING = "NA"


def _unbox(value: object) -> object:
    """Unwrap a numpy scalar or a length-one array, as R writes single values."""
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.item() if value.size == 1 else value
    if isinstance(value, str) or value is None:
        return value
    if isinstance(value, Sequence):
        items = cast(Sequence[object], value)
        return items[0] if len(items) == 1 else value
    return value


def _text(value: object, field: str) -> str:
    """Read a required string, however HDF5 handed it over."""
    unboxed = _unbox(value)
    if isinstance(unboxed, str):
        return unboxed
    if isinstance(unboxed, bytes):
        return unboxed.decode()
    raise ValueError(f"Expected a string for {field}, got {value!r}")


def _optional_text(value: object, field: str) -> str | None:
    """Read a role the backend may leave unset, such as a missing p-value."""
    if value is None:
        return None
    try:
        text = _text(value, field)
    except ValueError:
        return None
    if not text or text.lower() in {"na", "nan", "none"}:
        return None
    return text


def _names(value: object, field: str) -> list[str]:
    """Read an ordered list of column or key names."""
    if isinstance(value, str):
        return [value]
    if isinstance(value, np.ndarray):
        flattened = cast(list[object], value.reshape(-1).tolist())
        return [_text(item, field) for item in flattened]
    if isinstance(value, Sequence):
        return [_text(item, field) for item in cast(Sequence[object], value)]
    raise ValueError(f"Missing {field} in the prolfquapp AnnData metadata")


@dataclass(frozen=True, slots=True)
class ContrastRoles:
    """Which column carries which contrast fact, per the artifact itself."""

    subject_id: str
    contrast_col: str
    effect_col: str
    fdr_col: str
    model_name_col: str
    pvalue_col: str | None
    score_col: str | None
    directional: bool

    @classmethod
    def from_metadata(cls, metadata: Mapping[str, object]) -> ContrastRoles:
        configuration = metadata.get("contrast_configuration")
        if not isinstance(configuration, Mapping):
            raise ValueError(
                "AnnData uns['prolfquapp'] has no contrast_configuration: "
                "the column roles are missing, so the ranks cannot be resolved"
            )
        return cls(
            subject_id=_text(configuration.get("subject_id"), "subject_id"),
            contrast_col=_text(configuration.get("contrast_col"), "contrast_col"),
            effect_col=_text(configuration.get("effect_col"), "effect_col"),
            fdr_col=_text(configuration.get("fdr_col"), "fdr_col"),
            model_name_col=_text(configuration.get("model_name_col"), "model_name_col"),
            pvalue_col=_optional_text(configuration.get("pvalue_col"), "pvalue_col"),
            score_col=_optional_text(configuration.get("score_col"), "score_col"),
            directional=bool(_unbox(configuration.get("significance_directional"))),
        )

    @property
    def tests_differences(self) -> bool:
        """Whether the backend produced a p-value (prolfqua's `has_pvalue()`)."""
        return self.pvalue_col is not None

    def rank_score(self) -> pl.Expr:
        """The rank score, as prolfquapp's `.write_GSEA()` asks for it.

        The backend's test statistic, which is signed and unbounded. A backend
        that reports no p-value scores a bounded probability instead
        (SAINTexpress: SaintScore), which carries no direction and cannot order
        a ranked list, so its effect size is the rank.
        """
        if self.pvalue_col is None or self.score_col is None:
            return pl.col(self.effect_col).alias("score")
        return pl.col(self.score_col).alias("score")


@dataclass(frozen=True, slots=True)
class DeaArtifact:
    """Contrast rows and column roles read from one DEA result artifact."""

    roles: ContrastRoles
    rows: pl.DataFrame
    identifier: str

    @property
    def contrasts(self) -> list[str]:
        return self.rows.get_column(self.roles.contrast_col).unique().to_list()


def is_dea_artifact(path: Path) -> bool:
    """Whether this `.h5ad` is a DEA artifact whose column roles can be read.

    prolfquapp has written `AnnData.h5ad` beside its report since 2.9.0, but
    only 2.10.0 began recording `contrast_configuration`. Without those roles
    there is nothing to say which column is the contrast or the effect, so an
    older file is not a rank source and the archive's XLSX sheet is read
    instead. Reads the metadata alone, leaving the matrices on disk.
    """
    adata = ad.read_h5ad(path, backed="r")
    try:
        metadata = adata.uns.get(_PROLFQUAPP_UNS)
        if not isinstance(metadata, Mapping):
            return False
        has_roles = isinstance(metadata.get("contrast_configuration"), Mapping)
        return (
            has_roles
            and _optional_text(metadata.get("artifact_type"), "artifact_type") == _DEA_ARTIFACT_TYPE
        )
    finally:
        adata.file.close()


def read_dea_artifact(path: Path) -> DeaArtifact:
    """Read the contrast rows of a prolfquapp DEA AnnData file."""
    adata = ad.read_h5ad(path)
    metadata = adata.uns.get(_PROLFQUAPP_UNS)
    if not isinstance(metadata, Mapping):
        raise ValueError(
            f"AnnData was not written by prolfquapp (no uns['{_PROLFQUAPP_UNS}']): {path}"
        )
    artifact_type = _optional_text(metadata.get("artifact_type"), "artifact_type")
    if artifact_type != _DEA_ARTIFACT_TYPE:
        raise ValueError(
            f"AnnData is not a prolfquapp {_DEA_ARTIFACT_TYPE} artifact: got {artifact_type!r}"
        )

    roles = ContrastRoles.from_metadata(metadata)
    features = _feature_table(adata, roles)
    identifier = _text(metadata.get(IDENTIFIER_KEY), IDENTIFIER_KEY)
    if identifier not in features.columns:
        raise ValueError(
            f"AnnData {IDENTIFIER_KEY} names {identifier!r}, which is not a "
            f"feature annotation column: {features.columns}"
        )

    layout = _varm_layout(metadata)
    frames = [
        _with_features(features, layout.contrast_values(adata, key))
        for key in _names(metadata.get("varm_key_order"), "varm_key_order")
        if _is_contrast_key(key)
    ]
    if not frames:
        raise ValueError(f"AnnData carries no contrast results: {path}")
    return DeaArtifact(
        roles=roles,
        rows=pl.concat(frames, how="vertical_relaxed"),
        identifier=identifier,
    )


def _feature_table(adata: ad.AnnData, roles: ContrastRoles) -> pl.DataFrame:
    """The `var` annotation, keyed by the feature axis.

    Features whose annotation is missing carry R's "NA" text; restoring it to a
    real null lets them drop out of the ranks, as they do in prolfquapp.
    """
    # read_h5ad() builds an in-memory AnnData, whose var is a pandas frame.
    var = cast(pd.DataFrame, adata.var).reset_index(drop=False)
    frame = pl.from_pandas(var).with_columns(pl.col(pl.String).replace({_R_MISSING_STRING: None}))
    first = frame.columns[0]
    if roles.subject_id in frame.columns:
        return frame.drop(first)
    return frame.rename({first: roles.subject_id})


def _is_contrast_key(key: str) -> bool:
    """prolfquapp names a contrast frame with this prefix (its own spelling)."""
    return key.startswith("constrast_")


class _VarmLayout(Protocol):
    """How one artifact stores a contrast's columns in `varm`."""

    def contrast_values(self, adata: ad.AnnData, key: str) -> pl.DataFrame:
        """One contrast's columns, row-aligned with `var`."""
        ...


@dataclass(frozen=True, slots=True)
class _MatrixLayout:
    """prolfquapp 2.10.x: a numeric matrix, its column names kept in `uns`."""

    metadata: Mapping[str, object]

    def contrast_values(self, adata: ad.AnnData, key: str) -> pl.DataFrame:
        columns = _names(_group_entry(self.metadata, "varm_columns", key), f"varm_columns[{key}]")
        values = cast(Mapping[str, NDArray[np.float64]], adata.varm)[key]
        if values.shape[1] != len(columns):
            raise ValueError(
                f"AnnData varm[{key}] has {values.shape[1]} columns but "
                f"varm_columns names {len(columns)}"
            )
        frame = pl.DataFrame(
            {name: values[:, index] for index, name in enumerate(columns)},
        )
        annotations = _group_entry(self.metadata, "varm_annotations", key)
        if isinstance(annotations, Mapping):
            frame = frame.with_columns(
                [pl.Series(name, list(column)) for name, column in annotations.items()]
            )
        return frame


class _FrameLayout:
    """prolfquapp 2.11.0+: a data frame indexed by the feature ids."""

    def contrast_values(self, adata: ad.AnnData, key: str) -> pl.DataFrame:
        values = cast(Mapping[str, object], adata.varm)[key]
        if not isinstance(values, pd.DataFrame):
            raise ValueError(
                f"AnnData varm[{key}] is not a data frame, and no varm_columns names its columns"
            )
        return pl.from_pandas(values.reset_index(drop=True))


def _varm_layout(metadata: Mapping[str, object]) -> _VarmLayout:
    """The layout the writer used: only 2.10.x names the matrix columns in `uns`."""
    if isinstance(metadata.get("varm_columns"), Mapping):
        return _MatrixLayout(metadata)
    return _FrameLayout()


def _with_features(features: pl.DataFrame, values: pl.DataFrame) -> pl.DataFrame:
    """One contrast's rows, feature annotation joined on the feature axis.

    A 2.11.0 contrast frame repeats the feature keys `var` already carries;
    the annotation's copy is kept.
    """
    return pl.concat([features, values.select(pl.exclude(features.columns))], how="horizontal")


def _group_entry(metadata: Mapping[str, object], group: str, key: str) -> object:
    """Read one keyed entry from an `uns` group written per contrast."""
    entries = metadata.get(group)
    if not isinstance(entries, Mapping):
        return None
    return entries.get(key)
