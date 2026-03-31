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


@dataclass(frozen=True, slots=True)
class RankList:
    """Original ranked gene list submitted to STRING-DB for one contrast."""

    contrast: str
    entries: dict[str, float]  # input_label → score

    @property
    def n_genes(self) -> int:
        return len(self.entries)


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

    def mean_input_value(self, gene_pool: "GenePool") -> float:
        """Mean ranking score across gene hits — signed pseudo-NES."""
        hits = [gene_pool[gid] for gid in self.gene_ids if gid in gene_pool]
        if not hits:
            return 0.0
        return sum(h.input_value for h in hits) / len(hits)

    def rank_nes(self, gene_pool: "GenePool", n_input_genes: int) -> float:
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
    rank_lists: dict[str, RankList]

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
    These are protected during splitting, matching the logic in network.py.
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


def parse_gsea_tsv(
    path: Path,
    *,
    contrast: str | None = None,
    categories: set[str] | None = None,
) -> dict[str, CategoryGSEA]:
    """Parse a single-contrast STRING-DB GSEA TSV into CategoryGSEA objects.

    Builds one shared gene pool from all rows (all categories), then
    passes it to each category. All CategoryGSEA instances for the same
    contrast share the same pool by reference.

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

    # Build shared pool from ALL categories before filtering
    gene_pool = _build_gene_pool(df)

    if categories is not None:
        df = df.filter(pl.col("category").is_in(list(categories)))

    result: dict[str, CategoryGSEA] = {}
    for cat_name, cat_df in df.group_by("category"):
        name = cat_name[0]
        result[name] = _parse_category_group(cat_df, contrast, name, gene_pool)

    return result


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
