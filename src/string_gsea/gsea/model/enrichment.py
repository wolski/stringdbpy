# pyright: reportUnknownMemberType=false
"""Typed enrichment domain objects and STRING-DB TSV parsers."""

from __future__ import annotations

import json
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import cast

import polars as pl

from string_gsea.gsea.model.json_boundary import (
    json_array,
    json_object,
    required_integer,
    required_number,
    required_string,
    string_map,
)
from string_gsea.gsea.model.ranks import (
    GeneHit,
    GenePool,
    JsonObject,
    RankList,
    RankListCollection,
)


@dataclass(frozen=True, slots=True)
class TermGSEA:
    """One enriched term from STRING-DB GSEA.

    Gene detail is not embedded — terms hold protein IDs that
    reference the shared GenePool on the parent CategoryGSEA.

    STRING TSV column mapping:
        term_id          ← termID              (e.g. "GO:0006364")
        category         ← category            (e.g. "GO Process", "KEGG", "Reactome")
        description      ← termDescription     (e.g. "rRNA processing")
        enrichment_score ← enrichmentScore     (unsigned KS statistic)
        direction        ← direction           ("top", "bottom", "both ends")
        fdr              ← falseDiscoveryRate  (FDR-adjusted p-value)
        method           ← method              (e.g. "ks")
        genes_mapped     ← genesMapped         (input genes hitting this term)
        genes_in_set     ← genesInSet          (total genes annotated to this term)
        gene_ids         ← proteinIDs          (split from comma-separated string)
    """

    term_id: str
    category: str
    description: str
    enrichment_score: float
    direction: str
    fdr: float
    method: str
    genes_mapped: int
    genes_in_set: int
    gene_ids: tuple[str, ...]
    leading_edge_ids: tuple[str, ...] | None = None

    @property
    def gene_ratio(self) -> float:
        """Fraction of term genes found in input."""
        return self.genes_mapped / self.genes_in_set if self.genes_in_set > 0 else 0.0

    def mean_input_value(self, gene_pool: GenePool) -> float:
        """Mean ranking score across gene hits — signed pseudo-NES."""
        hits = [gene_pool[gid] for gid in self.gene_ids if gid in gene_pool]
        if not hits:
            return 0.0
        return sum(h.input_value for h in hits) / len(hits)

    def to_dict(self) -> JsonObject:
        data: JsonObject = {
            "term_id": self.term_id,
            "category": self.category,
            "description": self.description,
            "enrichment_score": self.enrichment_score,
            "direction": self.direction,
            "fdr": self.fdr,
            "method": self.method,
            "genes_mapped": self.genes_mapped,
            "genes_in_set": self.genes_in_set,
            "gene_ids": list(self.gene_ids),
        }
        if self.leading_edge_ids is not None:
            data["leading_edge_ids"] = list(self.leading_edge_ids)
        return data

    @classmethod
    def from_dict(cls, data: Mapping[str, object]) -> TermGSEA:
        return cls(
            term_id=required_string(data.get("term_id"), "term_id"),
            category=required_string(data.get("category"), "category"),
            description=required_string(data.get("description"), "description"),
            enrichment_score=required_number(data.get("enrichment_score"), "enrichment_score"),
            direction=required_string(data.get("direction"), "direction"),
            fdr=required_number(data.get("fdr"), "fdr"),
            method=required_string(data.get("method"), "method"),
            genes_mapped=required_integer(data.get("genes_mapped"), "genes_mapped"),
            genes_in_set=required_integer(data.get("genes_in_set"), "genes_in_set"),
            gene_ids=tuple(
                required_string(item, "gene_ids entry")
                for item in json_array(data.get("gene_ids"), "gene_ids")
            ),
            leading_edge_ids=(
                tuple(
                    required_string(item, "leading_edge_ids entry")
                    for item in json_array(data["leading_edge_ids"], "leading_edge_ids")
                )
                if "leading_edge_ids" in data
                else None
            ),
        )

    def rank_nes(self, gene_pool: GenePool, n_input_genes: int) -> float:
        """Ad-hoc NES based on mean rank position.

        Normalized to [-1, 1]: +1 = all genes at top of list,
        -1 = all at bottom, 0 = uniformly distributed.
        """
        hits = [gene_pool[gid] for gid in self.gene_ids if gid in gene_pool]
        if not hits:
            return 0.0
        mean_rank = sum(h.rank for h in hits) / len(hits)
        return 1.0 - 2.0 * mean_rank / (n_input_genes + 1)


@dataclass(frozen=True, slots=True)
class CategoryGSEA:
    """Enrichment results for one gene set category in one contrast.

    Holds a reference to the shared gene pool so that terms can
    resolve their gene detail without embedding copies.
    """

    category: str
    contrast: str
    gene_pool: GenePool
    terms: tuple[TermGSEA, ...]
    # Native R plotting payload, retained without interpreting S4-specific data.
    gsea_result: JsonObject | None = None

    def __post_init__(self) -> None:
        for term in self.terms:
            missing = [gid for gid in term.gene_ids if gid not in self.gene_pool]
            if missing:
                raise ValueError(f"Term {term.term_id}: unknown protein_ids: {missing[:3]}")

    def to_dict(self) -> JsonObject:
        """Serialize category (without gene_pool — stored at contrast level)."""
        data: JsonObject = {
            "category": self.category,
            "contrast": self.contrast,
            "terms": [t.to_dict() for t in self.terms],
        }
        if self.gsea_result is not None:
            data["gsea_result"] = self.gsea_result
        return data

    @classmethod
    def from_dict(cls, data: Mapping[str, object], gene_pool: GenePool) -> CategoryGSEA:
        return cls(
            category=required_string(data.get("category"), "category"),
            contrast=required_string(data.get("contrast"), "contrast"),
            gene_pool=gene_pool,
            terms=tuple(
                TermGSEA.from_dict(json_object(term, "term"))
                for term in json_array(data.get("terms"), "terms")
            ),
            gsea_result=(
                cast(JsonObject, json_object(data["gsea_result"], "gsea_result"))
                if "gsea_result" in data
                else None
            ),
        )


@dataclass(frozen=True, slots=True)
class MultiCategoryGSEA:
    """All categories for one contrast (= what STRING returns for one submission).

    The ``gene_pool`` is the full ranked gene list for the contrast and is
    shared by every child ``CategoryGSEA``. It is owned here (not derived from
    a category) so a contrast with no enriched terms still serializes cleanly.
    """

    contrast: str
    gene_pool: GenePool
    categories: dict[str, CategoryGSEA]

    def to_dict(self) -> JsonObject:
        # gene_pool is stored once per contrast (all categories share it)
        return {
            "contrast": self.contrast,
            "gene_pool": self.gene_pool.to_dict(),
            "categories": {name: cat.to_dict() for name, cat in self.categories.items()},
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, object]) -> MultiCategoryGSEA:
        gene_pool = GenePool.from_dict(json_object(data.get("gene_pool"), "gene_pool"))
        raw_categories = json_object(data.get("categories"), "categories")
        categories = {
            name: CategoryGSEA.from_dict(json_object(category, f"categories.{name}"), gene_pool)
            for name, category in raw_categories.items()
        }
        return cls(
            contrast=required_string(data.get("contrast"), "contrast"),
            gene_pool=gene_pool,
            categories=categories,
        )


@dataclass(frozen=True, slots=True)
class MultiContrastGSEA:
    """One category across multiple contrasts (= for cross-contrast plots)."""

    category: str
    contrasts: dict[str, CategoryGSEA]


@dataclass(frozen=True, slots=True)
class RunMetadata:
    """Parameters used for a STRING GSEA run.

    Captures everything needed to understand how results were produced,
    excluding sensitive fields (api_key).
    """

    workunit_id: str
    species: int
    fdr: float
    ge_enrichment_rank_direction: int
    caller_identity: str
    api_base_url: str

    def to_dict(self) -> JsonObject:
        return {
            "workunit_id": self.workunit_id,
            "species": self.species,
            "fdr": self.fdr,
            "ge_enrichment_rank_direction": self.ge_enrichment_rank_direction,
            "caller_identity": self.caller_identity,
            "api_base_url": self.api_base_url,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, object]) -> RunMetadata:
        return cls(
            workunit_id=required_string(data.get("workunit_id"), "workunit_id"),
            species=required_integer(data.get("species"), "species"),
            fdr=required_number(data.get("fdr"), "fdr"),
            ge_enrichment_rank_direction=required_integer(
                data.get("ge_enrichment_rank_direction"), "ge_enrichment_rank_direction"
            ),
            caller_identity=required_string(data.get("caller_identity"), "caller_identity"),
            api_base_url=required_string(data.get("api_base_url"), "api_base_url"),
        )


@dataclass(frozen=True, slots=True)
class GSEAResult:
    """Complete result: all contrasts x all categories.

    The data dict is keyed by contrast name. Use the slicing methods
    to get multi-category or multi-contrast views.
    """

    data: dict[str, MultiCategoryGSEA]
    rank_lists: dict[str, RankList]
    metadata: RunMetadata | None = None
    links: dict[str, str] | None = None  # contrast_key -> page_url

    @property
    def contrast_names(self) -> list[str]:
        return list(self.data.keys())

    @property
    def category_names(self) -> list[str]:
        cats: set[str] = set()
        for mc in self.data.values():
            cats.update(mc.categories.keys())
        return sorted(cats)

    def get_category(self, contrast: str, category: str) -> CategoryGSEA:
        return self.data[contrast].categories[category]

    def get_multi_category(self, contrast: str) -> MultiCategoryGSEA:
        return self.data[contrast]

    def get_multi_contrast(self, category: str) -> MultiContrastGSEA:
        contrasts = {
            name: mc.categories[category]
            for name, mc in self.data.items()
            if category in mc.categories
        }
        return MultiContrastGSEA(category=category, contrasts=contrasts)

    def to_dict(self) -> JsonObject:
        data: JsonObject = {
            "data": {name: mc.to_dict() for name, mc in self.data.items()},
            "rank_lists": {name: rl.to_dict() for name, rl in self.rank_lists.items()},
        }
        if self.metadata is not None:
            data["metadata"] = self.metadata.to_dict()
        if self.links is not None:
            data["links"] = dict(self.links)
        return data

    @classmethod
    def from_dict(cls, raw: Mapping[str, object]) -> GSEAResult:
        raw_data = json_object(raw.get("data"), "data")
        data = {
            name: MultiCategoryGSEA.from_dict(json_object(value, f"data.{name}"))
            for name, value in raw_data.items()
        }
        raw_rank_lists = json_object(raw.get("rank_lists"), "rank_lists")
        rank_lists = {
            name: RankList.from_dict(json_object(value, f"rank_lists.{name}"))
            for name, value in raw_rank_lists.items()
        }
        metadata = (
            RunMetadata.from_dict(json_object(raw["metadata"], "metadata"))
            if "metadata" in raw
            else None
        )
        raw_links = raw.get("links")
        links = None if raw_links is None else string_map(raw_links, "links")
        return cls(data=data, rank_lists=rank_lists, metadata=metadata, links=links)

    def to_json(self, path: Path) -> Path:
        """Serialize to JSON file."""
        with open(path, "w") as f:
            json.dump(self.to_dict(), f)
        return path

    @classmethod
    def from_json(cls, path: Path) -> GSEAResult:
        """Deserialize from JSON file."""
        with open(path) as f:
            payload: object = json.load(f)
        return cls.from_dict(json_object(payload, "GSEA result"))

    def to_polars_long(self) -> pl.DataFrame:
        """Convert to long-format Polars DataFrame matching STRING TSV schema.

        Produces one row per term (across all contrasts and categories) with
        the original STRING TSV columns reconstructed from the typed models,
        plus computed columns: directionNR, num_contrasts.
        """
        rows: list[dict[str, object]] = []
        for contrast_name, mc in self.data.items():
            for cat_name, cat in mc.categories.items():
                pool = cat.gene_pool
                for term in cat.terms:
                    # Reconstruct comma-separated protein columns from gene pool
                    pids = list(term.gene_ids)
                    hits = [pool.get(pid) for pid in pids]
                    rows.append(
                        {
                            "contrast": contrast_name,
                            "category": cat_name,
                            "termID": term.term_id,
                            "termDescription": term.description,
                            "enrichmentScore": term.enrichment_score,
                            "direction": term.direction,
                            "falseDiscoveryRate": term.fdr,
                            "method": term.method,
                            "genesMapped": term.genes_mapped,
                            "genesInSet": term.genes_in_set,
                            "proteinIDs": ",".join(pids),
                            "proteinLabels": ",".join(
                                h.label if h else pid for pid, h in zip(pids, hits, strict=True)
                            ),
                            "proteinInputLabels": ",".join(
                                h.input_label if h else pid
                                for pid, h in zip(pids, hits, strict=True)
                            ),
                            "proteinInputValues": ",".join(
                                str(h.input_value) if h else "" for h in hits
                            ),
                            "proteinRanks": ",".join(str(h.rank) if h else "" for h in hits),
                        }
                    )

        df = pl.DataFrame(rows)
        # Add directionNR
        df = df.with_columns(
            pl.when(pl.col("direction") == "top")
            .then(1)
            .when(pl.col("direction") == "bottom")
            .then(-1)
            .otherwise(0)
            .alias("directionNR")
        )
        # Add num_contrasts
        grouped = df.group_by(["category", "termID"]).agg(
            pl.count("contrast").alias("num_contrasts")
        )
        df = df.join(grouped, on=["category", "termID"], how="inner")
        return df

    def mapping_efficiency(self, contrast: str) -> float:
        """Fraction of submitted genes that STRING mapped (appeared in results).

        All categories for a contrast share the same gene pool, so we
        just need one category's pool to get the full set of mapped genes.
        """
        rank_list = self.rank_lists[contrast]
        mc = self.data[contrast]
        # All categories share the same pool — pick any one
        any_cat = next(iter(mc.categories.values()))
        n_mapped = len(any_cat.gene_pool)
        return n_mapped / rank_list.n_genes if rank_list.n_genes > 0 else 0.0


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------


def _split_protein_labels(raw: str) -> list[str]:
    """Split comma-separated proteinLabels, protecting embedded commas.

    Some protein labels contain commas followed by 1-2 digits (e.g. "HIST1H2BK,1").
    These are protected during splitting.
    """
    protected = re.sub(r",(\d{1,2},)", r"§COMMA§\1", raw)
    parts = protected.split(",")
    return [p.replace("§COMMA§", ",") for p in parts]


def _build_gene_pool(df: pl.DataFrame) -> GenePool:
    """Build a shared gene pool from all rows in a contrast's DataFrame."""
    entries: dict[str, GeneHit] = {}
    for row in df.iter_rows(named=True):
        protein_ids = row["proteinIDs"].split(",")
        labels = _split_protein_labels(row["proteinLabels"])
        input_labels = row["proteinInputLabels"].split(",")
        input_values = row["proteinInputValues"].split(",")
        ranks = row["proteinRanks"].split(",")

        for pid, label, ilabel, ival, rank in zip(
            protein_ids, labels, input_labels, input_values, ranks, strict=True
        ):
            if pid not in entries:
                entries[pid] = GeneHit(
                    protein_id=pid,
                    label=label.strip(),
                    input_label=ilabel.strip(),
                    input_value=float(ival),
                    rank=int(float(rank)),
                )
    return GenePool(entries=entries)


def _parse_category_group(
    df: pl.DataFrame, contrast: str, category: str, gene_pool: GenePool
) -> CategoryGSEA:
    """Build a CategoryGSEA from rows for one category, using a shared gene pool."""
    terms: list[TermGSEA] = []

    for row in df.iter_rows(named=True):
        protein_ids = row["proteinIDs"].split(",")
        terms.append(
            TermGSEA(
                term_id=row["termID"],
                category=category,
                description=row["termDescription"],
                enrichment_score=float(row["enrichmentScore"]),
                direction=row["direction"],
                fdr=float(row["falseDiscoveryRate"]),
                method=row["method"],
                genes_mapped=int(row["genesMapped"]),
                genes_in_set=int(row["genesInSet"]),
                gene_ids=tuple(protein_ids),
            )
        )

    return CategoryGSEA(
        category=category,
        contrast=contrast,
        gene_pool=gene_pool,
        terms=tuple(terms),
    )


def _parse_gsea_df(
    df: pl.DataFrame,
    contrast: str,
    categories: set[str] | None = None,
) -> tuple[dict[str, CategoryGSEA], GenePool]:
    """Core parsing logic: DataFrame → (dict of CategoryGSEA, shared gene pool).

    Builds one shared gene pool from all rows (all categories), then
    passes it to each category. The gene pool is returned alongside the
    categories so it can be owned at the contrast level — a contrast with
    no enriched terms yields an empty category dict but a valid gene pool.
    """
    gene_pool = _build_gene_pool(df)

    if categories is not None:
        df = df.filter(pl.col("category").is_in(list(categories)))

    result: dict[str, CategoryGSEA] = {}
    for cat_name, cat_df in df.group_by("category"):
        name = cat_name[0]
        result[name] = _parse_category_group(cat_df, contrast, name, gene_pool)

    return result, gene_pool


def parse_gsea_tsv(
    path: Path,
    *,
    contrast: str | None = None,
    categories: set[str] | None = None,
) -> tuple[dict[str, CategoryGSEA], GenePool]:
    """Parse a single-contrast STRING-DB GSEA TSV into CategoryGSEA objects.

    Args:
        path: Path to the TSV file (one contrast).
        contrast: Override contrast name. Default: the TSV filename.
        categories: If provided, only parse these categories. None means all.

    Returns:
        A ``(categories, gene_pool)`` tuple: the categories dict keyed by
        category name, and the shared gene pool for the contrast.
    """
    df = pl.read_csv(path, separator="\t")
    if contrast is None:
        contrast = Path(path).name
    return _parse_gsea_df(df, contrast, categories)


def parse_gsea_tsv_from_string(
    content: str,
    *,
    contrast: str,
    categories: set[str] | None = None,
) -> tuple[dict[str, CategoryGSEA], GenePool]:
    """Parse a single-contrast STRING-DB GSEA TSV from a string.

    Same as ``parse_gsea_tsv`` but reads from an in-memory string
    instead of a file path.

    Args:
        content: TSV content as a string.
        contrast: Contrast name for the parsed categories.
        categories: If provided, only parse these categories. None means all.

    Returns:
        A ``(categories, gene_pool)`` tuple: the categories dict keyed by
        category name, and the shared gene pool for the contrast.
    """
    import io

    df = pl.read_csv(io.StringIO(content), separator="\t")
    return _parse_gsea_df(df, contrast, categories)


def parse_rank_file(path: Path, *, contrast: str | None = None) -> RankList:
    """Parse a .rnk file (space/tab-separated: input_label score, no header)."""
    entries: dict[str, float] = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            entries[parts[0]] = float(parts[1])
    if contrast is None:
        contrast = path.stem
    return RankList(contrast=contrast, entries=entries)


def parse_gsea_tsv_dir(
    directory: Path,
    *,
    categories: set[str] | None = None,
) -> GSEAResult:
    """Parse all *_results.tsv files in a directory into a GSEAResult.

    Args:
        directory: Directory containing one TSV per contrast.
        categories: If provided, only parse these categories. None means all.

    Returns:
        GSEAResult with one MultiCategoryGSEA per contrast.
    """
    tsv_files = sorted(directory.glob("*_results.tsv"))
    if not tsv_files:
        raise FileNotFoundError(f"No *_results.tsv files found in {directory}")

    # Match rank files to contrasts by name prefix
    rnk_by_stem = {p.stem: p for p in directory.glob("*.rnk")}

    data: dict[str, MultiCategoryGSEA] = {}
    rank_lists: dict[str, RankList] = {}
    for tsv_path in tsv_files:
        contrast = tsv_path.name
        cat_dict, gene_pool = parse_gsea_tsv(tsv_path, contrast=contrast, categories=categories)
        data[contrast] = MultiCategoryGSEA(
            contrast=contrast, gene_pool=gene_pool, categories=cat_dict
        )

        # e.g. "Bait_NCP_pUbT12_results.tsv" → stem "Bait_NCP_pUbT12_results" → strip "_results"
        rnk_stem = tsv_path.stem.removesuffix("_results")
        if rnk_stem in rnk_by_stem:
            rank_lists[contrast] = parse_rank_file(rnk_by_stem[rnk_stem], contrast=contrast)

    return GSEAResult(data=data, rank_lists=rank_lists)


def parse_gsea_results(
    rank_lists: RankListCollection,
    tsv_content: dict[tuple[str, str], str],
    *,
    metadata: RunMetadata | None = None,
    links: dict[str, str] | None = None,
    categories: set[str] | None = None,
) -> GSEAResult:
    """Build a GSEAResult from in-memory downloaded TSV content.

    Args:
        rank_lists: The submitted rank lists (for mapping efficiency).
        tsv_content: Dict mapping ``(analysis, contrast)`` to TSV text,
            as downloaded by the injected GSEA gateway.
        metadata: Optional run parameters (config, workunit_id, species).
        categories: If provided, only parse these categories.

    Returns:
        GSEAResult with one MultiCategoryGSEA per contrast.
    """
    data: dict[str, MultiCategoryGSEA] = {}
    rl_dict: dict[str, RankList] = {}

    for (_analysis, contrast), tsv_text in tsv_content.items():
        contrast_key = f"{contrast}_results.tsv"
        cat_dict, gene_pool = parse_gsea_tsv_from_string(
            tsv_text, contrast=contrast_key, categories=categories
        )
        data[contrast_key] = MultiCategoryGSEA(
            contrast=contrast_key, gene_pool=gene_pool, categories=cat_dict
        )

        if contrast in rank_lists:
            rl = rank_lists[contrast]
            rl_dict[contrast_key] = RankList(contrast=contrast_key, entries=rl.entries)

    return GSEAResult(data=data, rank_lists=rl_dict, metadata=metadata, links=links)
