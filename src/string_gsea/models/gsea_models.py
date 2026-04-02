"""Typed domain objects for STRING-DB GSEA results and a TSV parser."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import polars as pl

if TYPE_CHECKING:
    from string_gsea.gsea_config import GSEAConfig

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

    def to_dict(self) -> dict:
        return {
            "protein_id": self.protein_id,
            "label": self.label,
            "input_label": self.input_label,
            "input_value": self.input_value,
            "rank": self.rank,
        }

    @classmethod
    def from_dict(cls, d: dict) -> GeneHit:
        return cls(**d)


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

    def __iter__(self):
        return iter(self.entries)

    def values(self):
        return self.entries.values()

    def get(self, protein_id: str, default=None):
        return self.entries.get(protein_id, default)

    def to_dict(self) -> dict:
        return {pid: hit.to_dict() for pid, hit in self.entries.items()}

    @classmethod
    def from_dict(cls, d: dict) -> GenePool:
        return cls(entries={pid: GeneHit.from_dict(h) for pid, h in d.items()})


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

    def to_dict(self) -> dict:
        return {"contrast": self.contrast, "entries": self.entries}

    @classmethod
    def from_dict(cls, d: dict) -> RankList:
        return cls(contrast=d["contrast"], entries=d["entries"])

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

    def __init__(self, analysis: str, rank_lists: list[RankList]):
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

    def __iter__(self):
        return iter(self._data.values())

    def items(self):
        """Yields (contrast_name, RankList) pairs."""
        return self._data.items()

    def values(self):
        return self._data.values()

    def keys(self):
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

    def to_dict(self) -> dict:
        return {
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

    @classmethod
    def from_dict(cls, d: dict) -> TermGSEA:
        return cls(**{**d, "gene_ids": tuple(d["gene_ids"])})

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

    def __post_init__(self):
        for term in self.terms:
            missing = [gid for gid in term.gene_ids if gid not in self.gene_pool]
            if missing:
                raise ValueError(f"Term {term.term_id}: unknown protein_ids: {missing[:3]}")

    def to_dict(self) -> dict:
        """Serialize category (without gene_pool — stored at contrast level)."""
        return {
            "category": self.category,
            "contrast": self.contrast,
            "terms": [t.to_dict() for t in self.terms],
        }

    @classmethod
    def from_dict(cls, d: dict, gene_pool: GenePool) -> CategoryGSEA:
        return cls(
            category=d["category"],
            contrast=d["contrast"],
            gene_pool=gene_pool,
            terms=tuple(TermGSEA.from_dict(t) for t in d["terms"]),
        )


@dataclass(frozen=True, slots=True)
class MultiCategoryGSEA:
    """All categories for one contrast (= what STRING returns for one submission)."""

    contrast: str
    categories: dict[str, CategoryGSEA]

    def to_dict(self) -> dict:
        # Store gene_pool once per contrast (all categories share it)
        any_cat = next(iter(self.categories.values()))
        return {
            "contrast": self.contrast,
            "gene_pool": any_cat.gene_pool.to_dict(),
            "categories": {name: cat.to_dict() for name, cat in self.categories.items()},
        }

    @classmethod
    def from_dict(cls, d: dict) -> MultiCategoryGSEA:
        gene_pool = GenePool.from_dict(d["gene_pool"])
        categories = {name: CategoryGSEA.from_dict(cat_d, gene_pool) for name, cat_d in d["categories"].items()}
        return cls(contrast=d["contrast"], categories=categories)


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

    def to_dict(self) -> dict:
        return {
            "workunit_id": self.workunit_id,
            "species": self.species,
            "fdr": self.fdr,
            "ge_enrichment_rank_direction": self.ge_enrichment_rank_direction,
            "caller_identity": self.caller_identity,
            "api_base_url": self.api_base_url,
        }

    @classmethod
    def from_dict(cls, d: dict) -> RunMetadata:
        return cls(
            workunit_id=d["workunit_id"],
            species=d["species"],
            fdr=d["fdr"],
            ge_enrichment_rank_direction=d["ge_enrichment_rank_direction"],
            caller_identity=d["caller_identity"],
            api_base_url=d["api_base_url"],
        )

    @classmethod
    def from_config(cls, config: GSEAConfig, *, workunit_id: str, species: int) -> RunMetadata:
        """Build from a GSEAConfig, excluding the api_key."""
        return cls(
            workunit_id=workunit_id,
            species=species,
            fdr=config.fdr,
            ge_enrichment_rank_direction=config.ge_enrichment_rank_direction,
            caller_identity=config.caller_identity,
            api_base_url=config.api_base_url,
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
        contrasts = {name: mc.categories[category] for name, mc in self.data.items() if category in mc.categories}
        return MultiContrastGSEA(category=category, contrasts=contrasts)

    def to_dict(self) -> dict:
        d: dict = {
            "data": {name: mc.to_dict() for name, mc in self.data.items()},
            "rank_lists": {name: rl.to_dict() for name, rl in self.rank_lists.items()},
        }
        if self.metadata is not None:
            d["metadata"] = self.metadata.to_dict()
        if self.links is not None:
            d["links"] = self.links
        return d

    @classmethod
    def from_dict(cls, d: dict) -> GSEAResult:
        data = {name: MultiCategoryGSEA.from_dict(mc_d) for name, mc_d in d["data"].items()}
        rank_lists = {name: RankList.from_dict(rl_d) for name, rl_d in d["rank_lists"].items()}
        metadata = RunMetadata.from_dict(d["metadata"]) if "metadata" in d else None
        links = d.get("links")
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
            return cls.from_dict(json.load(f))

    def to_polars_long(self) -> pl.DataFrame:
        """Convert to long-format Polars DataFrame matching STRING TSV schema.

        Produces one row per term (across all contrasts and categories) with
        the original STRING TSV columns reconstructed from the typed models,
        plus computed columns: directionNR, num_contrasts.
        """
        rows: list[dict] = []
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
                            "proteinLabels": ",".join(h.label if h else pid for pid, h in zip(pids, hits, strict=True)),
                            "proteinInputLabels": ",".join(
                                h.input_label if h else pid for pid, h in zip(pids, hits, strict=True)
                            ),
                            "proteinInputValues": ",".join(str(h.input_value) if h else "" for h in hits),
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
        grouped = df.group_by(["category", "termID"]).agg(pl.count("contrast").alias("num_contrasts"))
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

        for pid, label, ilabel, ival, rank in zip(protein_ids, labels, input_labels, input_values, ranks, strict=True):
            if pid not in entries:
                entries[pid] = GeneHit(
                    protein_id=pid,
                    label=label.strip(),
                    input_label=ilabel.strip(),
                    input_value=float(ival),
                    rank=int(float(rank)),
                )
    return GenePool(entries=entries)


def _parse_category_group(df: pl.DataFrame, contrast: str, category: str, gene_pool: GenePool) -> CategoryGSEA:
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
) -> dict[str, CategoryGSEA]:
    """Core parsing logic: DataFrame → dict of CategoryGSEA.

    Builds one shared gene pool from all rows (all categories), then
    passes it to each category.
    """
    gene_pool = _build_gene_pool(df)

    if categories is not None:
        df = df.filter(pl.col("category").is_in(list(categories)))

    result: dict[str, CategoryGSEA] = {}
    for cat_name, cat_df in df.group_by("category"):
        name = cat_name[0]
        result[name] = _parse_category_group(cat_df, contrast, name, gene_pool)

    return result


def parse_gsea_tsv(
    path: Path,
    *,
    contrast: str | None = None,
    categories: set[str] | None = None,
) -> dict[str, CategoryGSEA]:
    """Parse a single-contrast STRING-DB GSEA TSV into CategoryGSEA objects.

    Args:
        path: Path to the TSV file (one contrast).
        contrast: Override contrast name. Default: the TSV filename.
        categories: If provided, only parse these categories. None means all.

    Returns:
        Dict keyed by category name.
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
) -> dict[str, CategoryGSEA]:
    """Parse a single-contrast STRING-DB GSEA TSV from a string.

    Same as ``parse_gsea_tsv`` but reads from an in-memory string
    instead of a file path.

    Args:
        content: TSV content as a string.
        contrast: Contrast name for the parsed categories.
        categories: If provided, only parse these categories. None means all.

    Returns:
        Dict keyed by category name.
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
        cat_dict = parse_gsea_tsv(tsv_path, contrast=contrast, categories=categories)
        data[contrast] = MultiCategoryGSEA(contrast=contrast, categories=cat_dict)

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
            as populated by ``StringGSEAResults.download()``.
        metadata: Optional run parameters (config, workunit_id, species).
        categories: If provided, only parse these categories.

    Returns:
        GSEAResult with one MultiCategoryGSEA per contrast.
    """
    data: dict[str, MultiCategoryGSEA] = {}
    rl_dict: dict[str, RankList] = {}

    for (_analysis, contrast), tsv_text in tsv_content.items():
        contrast_key = f"{contrast}_results.tsv"
        cat_dict = parse_gsea_tsv_from_string(tsv_text, contrast=contrast_key, categories=categories)
        data[contrast_key] = MultiCategoryGSEA(contrast=contrast_key, categories=cat_dict)

        if contrast in rank_lists:
            rl = rank_lists[contrast]
            rl_dict[contrast_key] = RankList(contrast=contrast_key, entries=rl.entries)

    return GSEAResult(data=data, rank_lists=rl_dict, metadata=metadata, links=links)
