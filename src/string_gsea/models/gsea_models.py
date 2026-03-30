"""Typed domain objects for STRING-DB GSEA results and a TSV parser."""

import re
from dataclasses import dataclass
from pathlib import Path

import polars as pl

# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class GeneHit:
    """A single gene/protein from the input ranking.

    Stored once per contrast in a shared pool (GenePool).
    Terms reference genes by protein_id.
    """

    protein_id: str
    label: str
    input_label: str
    input_value: float
    rank: int


GenePool = dict[str, GeneHit]


@dataclass(frozen=True, slots=True)
class TermGSEA:
    """One enriched term from STRING-DB GSEA.

    Gene detail is not embedded — terms hold protein IDs that
    reference the shared GenePool on the parent CategoryGSEA.
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


@dataclass(frozen=True, slots=True)
class MultiCategoryGSEA:
    """All categories for one contrast (= what STRING returns for one submission)."""

    contrast: str
    categories: dict[str, CategoryGSEA]


@dataclass(frozen=True, slots=True)
class MultiContrastGSEA:
    """One category across multiple contrasts (= for cross-contrast plots)."""

    category: str
    contrasts: dict[str, CategoryGSEA]


@dataclass(frozen=True, slots=True)
class GSEAResult:
    """Complete result: all contrasts x all categories.

    The data dict is keyed by contrast name. Use the slicing methods
    to get multi-category or multi-contrast views.
    """

    data: dict[str, MultiCategoryGSEA]

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


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------


def _split_protein_labels(raw: str) -> list[str]:
    """Split comma-separated proteinLabels, protecting embedded commas.

    Some protein labels contain commas followed by 1-2 digits (e.g. "HIST1H2BK,1").
    These are protected during splitting, matching the logic in network.py.
    """
    protected = re.sub(r",(\d{1,2},)", r"§COMMA§\1", raw)
    parts = protected.split(",")
    return [p.replace("§COMMA§", ",") for p in parts]


def _parse_category_group(df: pl.DataFrame, contrast: str, category: str) -> CategoryGSEA:
    """Build a CategoryGSEA from a DataFrame containing rows for one category."""
    gene_pool: GenePool = {}
    terms: list[TermGSEA] = []

    for row in df.iter_rows(named=True):
        protein_ids = row["proteinIDs"].split(",")
        labels = _split_protein_labels(row["proteinLabels"])
        input_labels = row["proteinInputLabels"].split(",")
        input_values = row["proteinInputValues"].split(",")
        ranks = row["proteinRanks"].split(",")

        for pid, label, ilabel, ival, rank in zip(protein_ids, labels, input_labels, input_values, ranks, strict=True):
            if pid not in gene_pool:
                gene_pool[pid] = GeneHit(
                    protein_id=pid,
                    label=label.strip(),
                    input_label=ilabel.strip(),
                    input_value=float(ival),
                    rank=int(float(rank)),
                )

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

    if categories is not None:
        df = df.filter(pl.col("category").is_in(list(categories)))

    result: dict[str, CategoryGSEA] = {}
    for cat_name, cat_df in df.group_by("category"):
        name = cat_name[0]
        result[name] = _parse_category_group(cat_df, contrast, name)

    return result


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

    data: dict[str, MultiCategoryGSEA] = {}
    for tsv_path in tsv_files:
        contrast = tsv_path.name
        cat_dict = parse_gsea_tsv(tsv_path, contrast=contrast, categories=categories)
        data[contrast] = MultiCategoryGSEA(contrast=contrast, categories=cat_dict)

    return GSEAResult(data=data)
