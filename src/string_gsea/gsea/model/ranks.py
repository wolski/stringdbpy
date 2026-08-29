"""Typed ranked-gene domain objects for STRING-GSEA."""

from __future__ import annotations

from collections.abc import ItemsView, Iterator, KeysView, Mapping, ValuesView
from dataclasses import dataclass

import polars as pl

from string_gsea.gsea.model.json_boundary import (
    json_object,
    required_integer,
    required_number,
    required_string,
)

type JsonValue = str | int | float | bool | list[JsonValue] | dict[str, JsonValue] | None
type JsonObject = dict[str, JsonValue]

# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class GeneHit:
    """A single gene/protein from the input ranking.

    Stored once per contrast in a shared pool (GenePool).
    Terms reference genes by protein_id.

    STRING TSV column mapping:
        protein_id   ← proteinIDs          (STRING internal ID, e.g. "9606.ENSP00000007390")
        label        ← proteinLabels       (gene symbol, e.g. "TSR3")
        input_label  ← proteinInputLabels  (original submitted ID, e.g. "Q9UJK0")
        input_value  ← proteinInputValues  (ranking score, e.g. fold change or t-statistic)
        rank         ← proteinRanks        (rank position in the input list)
    """

    protein_id: str
    label: str
    input_label: str
    input_value: float
    rank: int

    def to_dict(self) -> JsonObject:
        return {
            "protein_id": self.protein_id,
            "label": self.label,
            "input_label": self.input_label,
            "input_value": self.input_value,
            "rank": self.rank,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, object]) -> GeneHit:
        return cls(
            protein_id=required_string(data.get("protein_id"), "protein_id"),
            label=required_string(data.get("label"), "label"),
            input_label=required_string(data.get("input_label"), "input_label"),
            input_value=required_number(data.get("input_value"), "input_value"),
            rank=required_integer(data.get("rank"), "rank"),
        )


@dataclass(frozen=True, slots=True)
class GenePool:
    """Shared pool of mapped genes for one contrast.

    Built once from all categories in a contrast's TSV. All CategoryGSEA
    instances for the same contrast share the same GenePool by reference.
    """

    entries: dict[str, GeneHit]  # keyed by protein_id

    @property
    def n_genes(self) -> int:
        return len(self.entries)

    def __contains__(self, protein_id: str) -> bool:
        return protein_id in self.entries

    def __getitem__(self, protein_id: str) -> GeneHit:
        return self.entries[protein_id]

    def __len__(self) -> int:
        return len(self.entries)

    def __iter__(self) -> Iterator[str]:
        return iter(self.entries)

    def values(self) -> ValuesView[GeneHit]:
        return self.entries.values()

    def get(self, protein_id: str, default: GeneHit | None = None) -> GeneHit | None:
        return self.entries.get(protein_id, default)

    def to_dict(self) -> JsonObject:
        return {pid: hit.to_dict() for pid, hit in self.entries.items()}

    @classmethod
    def from_dict(cls, data: Mapping[str, object]) -> GenePool:
        return cls(
            entries={
                protein_id: GeneHit.from_dict(json_object(hit, f"gene_pool.{protein_id}"))
                for protein_id, hit in data.items()
            }
        )


@dataclass(frozen=True, slots=True)
class RankList:
    """One ranked gene list for a single biological contrast.

    - ``contrast``: the biological comparison, e.g. ``"Treatment_vs_Control"``
    - ``entries``: mapping of input identifier → ranking score
    """

    contrast: str
    entries: dict[str, float]  # input_label → score

    @property
    def n_genes(self) -> int:
        return len(self.entries)

    def to_rnk_string(self) -> str:
        """Tab-separated rank string (no header), for STRING API and .rnk files."""
        return "\n".join(f"{label}\t{score}" for label, score in self.entries.items()) + "\n"

    def sample_identifiers(self, nr: int = 10) -> list[str]:
        """Return up to *nr* randomly sampled identifier keys."""
        import random

        keys = list(self.entries.keys())
        return random.sample(keys, min(nr, len(keys)))

    def to_dict(self) -> JsonObject:
        entries: JsonObject = dict(self.entries)
        return {"contrast": self.contrast, "entries": entries}

    @classmethod
    def from_dict(cls, data: Mapping[str, object]) -> RankList:
        raw_entries = json_object(data.get("entries"), "entries")
        return cls(
            contrast=required_string(data.get("contrast"), "contrast"),
            entries={
                identifier: required_number(value, f"entries.{identifier}")
                for identifier, value in raw_entries.items()
            },
        )

    @classmethod
    def from_polars(cls, df: pl.DataFrame, *, contrast: str) -> RankList:
        """Build from a 2-column DataFrame (identifier, score)."""
        entries = dict(
            zip(
                df.get_column(df.columns[0]).to_list(),
                df.get_column(df.columns[1]).cast(pl.Float64).to_list(),
                strict=True,
            )
        )
        return cls(contrast=contrast, entries=entries)


class RankListCollection:
    """Ranked gene lists for one analysis type across multiple contrasts.

    - ``analysis``: the filtering strategy that produced these lists, e.g.
      ``"pep_2_no_imputed"`` (≥2 peptides, no imputed), ``"pep_1"`` (all
      peptides), or ``"from_rnk"`` (raw .rnk file input).
    - Each ``RankList`` inside represents one biological contrast.

    Supports dict-like lookup by contrast name and iteration over
    ``RankList`` values.
    """

    def __init__(self, analysis: str, rank_lists: list[RankList]) -> None:
        self.analysis = analysis
        self._data: dict[str, RankList] = {}
        for rl in rank_lists:
            self._data[rl.contrast] = rl

    def add(self, rl: RankList) -> None:
        self._data[rl.contrast] = rl

    # --- dict-like interface ---------------------------------------------------

    def __getitem__(self, contrast: str) -> RankList:
        return self._data[contrast]

    def __contains__(self, contrast: str) -> bool:
        return contrast in self._data

    def __len__(self) -> int:
        return len(self._data)

    def __iter__(self) -> Iterator[RankList]:
        return iter(self._data.values())

    def items(self) -> ItemsView[str, RankList]:
        """Yields (contrast_name, RankList) pairs."""
        return self._data.items()

    def values(self) -> ValuesView[RankList]:
        return self._data.values()

    def keys(self) -> KeysView[str]:
        """Yields contrast names."""
        return self._data.keys()

    # --- accessors ------------------------------------------------------------

    @property
    def contrasts(self) -> list[str]:
        """All contrast names."""
        return sorted(self._data.keys())

    def first(self) -> RankList:
        """Return the first RankList (useful for species-detection fallback)."""
        return next(iter(self._data.values()))
