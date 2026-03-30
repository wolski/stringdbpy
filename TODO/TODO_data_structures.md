# TODO: Data Structures for Enrichment Results

Date: 2026-03-30

## Motivation

The pipeline currently passes untyped Polars DataFrames between stages, relying on implicit column conventions. This leads to:

1. **No documentation of what columns exist at each stage** — you have to read the code
2. **No validation** — misspelled column names or missing columns fail silently or with cryptic Polars errors
3. **Serialization scattered across files** — XLSX writing in `gsea_result_processor.py`, TSV/JSON in `string_ora_run.py`, no shared abstraction
4. **Difficulty comparing with other tools** — no clear mapping between STRING-DB's schema and the community's naming conventions

---

## Landscape: How Other Tools Model Enrichment Results

### Python packages

| Package | Result type | Typed? | Notes |
|---------|-------------|--------|-------|
| **GSEApy** | `pd.DataFrame` via `obj.res2d` | No | Convention columns: `Term`, `ES`, `NES`, `FDR q-val`, `Lead_genes` |
| **blitzgsea** | `pd.DataFrame` | No | Columns: `es`, `nes`, `pval`, `fdr`, `geneset_size`, `leading_edge` |
| **decoupler** | AnnData `.obsm` matrices | No (matrix slots) | Scores + pvals as obs x source matrices; no long-format schema |
| **goatools** | `GOEnrichmentRecord` (plain class) | Yes (attributes) | `GO`, `NS`, `enrichment`, `p_uncorrected`, `p_fdr_bh`, `study_items` |
| **g:Profiler** | `pd.DataFrame` | No | Notable: has `source` column (= category), `precision`, `recall` |

### R packages (reference)

| Package | Result type | Typed? | Notes |
|---------|-------------|--------|-------|
| **clusterProfiler** | S4 `enrichResult` / `gseaResult` | Yes | De facto standard: `ID`, `Description`, `NES`, `p.adjust`, `core_enrichment` |
| **fgsea** | `data.table` | No | `pathway`, `pval`, `padj`, `ES`, `NES`, `size`, `leadingEdge` |

### Key observations

1. **No standard interchange format exists** for enrichment results (unlike GMT for gene sets).
2. Most Python tools use **untyped DataFrames with convention columns**. Only goatools has a typed result object.
3. R's clusterProfiler S4 classes are the closest thing to a community standard schema.
4. **g:Profiler** is the best analog to STRING-DB's multi-category approach — it returns a `source` column (GO:BP, KEGG, REAC, etc.) comparable to STRING's `category`.

### What makes STRING-DB GSEA different

| Feature | STRING-DB | clusterProfiler / GSEApy | g:Profiler |
|---------|-----------|--------------------------|------------|
| Multi-category results | Yes (15 categories in one response) | No (one library per run) | Yes (`source` column) |
| Direction field | Explicit `direction` column: "top", "bottom", "both ends" | Encoded in sign of ES/NES | N/A (ORA only) |
| Enrichment score | `enrichmentScore` (unsigned KS statistic) | `ES` (signed) + `NES` (normalized) | N/A |
| Per-protein values | `proteinInputValues`, `proteinRanks` per term | `leading_edge` (gene list only) | `intersection` (gene list only) |
| FDR only | `falseDiscoveryRate` (no raw p-value) | Both raw and adjusted | Adjusted only |

The STRING API returns **unsigned enrichment scores** with a separate `direction` field, unlike classic GSEA which uses signed ES/NES. STRING also returns per-protein input values and ranks within each term — richer than leading-edge gene lists.

---

## Current Column Schema by Pipeline Stage

### Stage 1: Raw TSV from STRING API

```
category, termID, termDescription, genesMapped, genesInSet,
enrichmentScore, direction, falseDiscoveryRate, method,
proteinIDs, proteinLabels, proteinInputLabels,
proteinInputValues, proteinRanks
```

All Utf8 in TSV. Protein columns are comma-separated strings.

### Stage 2: Long XLSX (`_results_to_dataframe()`)

Adds:
- `contrast` (Utf8) — from filename
- `directionNR` (Int64) — 1=top, -1=bottom, 0=other
- `num_contrasts` (Int64) — count per (category, termID)


### Stage 3: After `filter_by_FDR()`

Rows filtered: `falseDiscoveryRate < threshold AND genesMapped > threshold`

### Stage 4: After `select_top_terms()`

Top N terms per (contrast, category) by FDR. No column changes.

### Stage 5: After `add_gene_ratio()`

Adds: `geneRatio` (Float64) = genesMapped / genesInSet

### Stage 6: After `explode_protein_columns()`

Comma-separated protein columns split to individual rows. `proteinInputValues` and `proteinRanks` cast to Float64.

### Stage 7: After `summarize_terms()`

Adds: `meanInputValues` (Float64) — mean proteinInputValues per (contrast, termID)

### Note on pseudo-NES

STRING returns unsigned `enrichmentScore` + a `direction` label. But we have `proteinInputValues` (the original ranking scores, e.g., fold changes or t-statistics) per term. We can compute a **mean input value** per term — this gives a signed score that naturally captures direction:

- Positive mean → genes in this term tend to be upregulated
- Negative mean → genes tend to be downregulated

This is already partially implemented as `meanInputValues` in `summarize_terms()`. It could serve as a pseudo-NES for downstream use (dotplots, sorting, coloring).

### ORA results (separate pipeline)

COMMENT: No need to adress this one at this point.

From `/json/enrichment` endpoint: `category`, `term`, `description`, `number_of_genes`, `number_of_genes_in_background`, `pvalue`, `fdr`, `inputGenes`

Different column names than GSEA — no shared schema.

---

## Naming Discussion: Container Hierarchy

The class hierarchy maps to how the data is naturally sliced:

| Level | What it holds | Class name |
|-------|--------------|------------|
| One term, one contrast | Single enrichment result | `TermGSEA` |
| One category, one contrast | All terms for e.g. "GO Process" in contrast A | `CategoryGSEA` |
| All categories, one contrast | Everything STRING returns for one submission | `MultiCategoryGSEA` |
| One category, all contrasts | e.g. "GO Process" across contrasts A, B, C | `MultiContrastGSEA` |
| Everything | All categories × all contrasts | `GSEAResult` |

**Naming convention:** Names describe the *shape* of the container (what axes it spans), not a domain concept. This was debated — `ContrastGSEA` (naming by what it *represents*: one contrast's full result) vs `MultiCategoryGSEA` (naming by *structure*: multiple categories). Decision: keep structural names for consistency across the hierarchy.

**Slicing model:** `GSEAResult` is the top-level container. The multi-* classes are *views/slices*, not separate storage:

- `GSEAResult` → pick one contrast → `MultiCategoryGSEA` (all categories for that contrast = what STRING returns)
- `GSEAResult` → pick one category → `MultiContrastGSEA` (that category across all contrasts = for cross-contrast dotplots)
- Either of those → pick one (category, contrast) → `CategoryGSEA`

---

## Proposed Architecture: Domain Object Model

### Design philosophy

Instead of passing flat DataFrames through the pipeline, decompose the long table into **typed domain objects** that map directly to visualization and analysis use cases. Inspired by clusterProfiler's `enrichResult` / `gseaResult` S4 classes, but adapted for STRING-DB's multi-category, multi-contrast structure.

The key insight: **each visualization function should accept exactly the data structure it needs.** A per-category ridge plot takes a `CategoryGSEA`. A cross-contrast dotplot takes a `MultiContrastGSEA`. The type signature documents what data the function expects.

### Module: `src/string_gsea/models/`

```
src/string_gsea/models/
├── __init__.py              # Public re-exports
├── gene_hit.py              # GeneHit — per-gene detail within a term
├── term_gsea.py             # TermGSEA — a single enriched term
├── category_gsea.py         # CategoryGSEA — one category, one contrast
├── containers.py            # MultiCategoryGSEA, MultiContrastGSEA
├── ora_result.py            # ORA equivalents
└── serializers.py           # XLSX, TSV, JSON, clusterProfiler-format export
```

---

## Phase 1: Core Pydantic Models — Bottom-Up

### 1.1 `GeneHit` — per-gene detail (shared pool)

A gene that appears in 5 significant terms has the same `input_value`, `rank`, and `label` in all 5 — this information comes from the input rank file and is global to the contrast, not specific to any term. Therefore, `GeneHit` objects live in a **shared gene pool** at the contrast level, and terms reference them by protein ID.

```python
from pydantic import BaseModel, Field


class GeneHit(BaseModel):
    """A single gene/protein from the input ranking.

    Stored once per contrast in a shared pool (GenePool).
    Terms reference genes by protein_id, not by embedding copies.
    """

    protein_id: str          # STRING ID (e.g., "9606.ENSP00000269305")
    label: str               # Gene symbol (e.g., "TP53")
    input_label: str         # Original input identifier
    input_value: float       # Ranking score (fold change, t-statistic, etc.)
    rank: float              # Rank position in the input list


# Type alias for the shared pool
GenePool = dict[str, GeneHit]   # keyed by protein_id
```

Currently these are comma-separated strings in the API response. `explode_protein_columns()` splits them into rows, exploding the DataFrame from ~2,883 to ~161,000 rows. With a shared gene pool, each gene is parsed **once** and referenced by ID from every term that contains it.

### 1.2 `TermGSEA` — a single enriched term (references genes by ID)

```python
class TermGSEA(BaseModel):
    """One enriched term from STRING-DB GSEA.

    Gene detail is not embedded — terms hold protein IDs that
    reference the shared GenePool on the parent container.
    """

    # Term identification
    term_id: str
    term_description: str

    # Enrichment statistics (from STRING API)
    enrichment_score: float       # Unsigned KS statistic
    direction: str                # "top", "bottom", "both ends"
    false_discovery_rate: float
    method: str                   # "ks" typically

    # Gene set sizes
    genes_mapped: int             # How many input genes hit this term
    genes_in_set: int             # Total genes annotated to this term

    # References into the shared gene pool (not embedded copies)
    gene_ids: list[str]           # protein_id values → look up in GenePool

    # Computed properties (that don't need the gene pool)
    @property
    def gene_ratio(self) -> float:
        """Fraction of term genes found in input."""
        return self.genes_mapped / self.genes_in_set if self.genes_in_set > 0 else 0.0

    @property
    def direction_nr(self) -> int:
        """Numeric direction: 1=top, -1=bottom, 0=both/other."""
        return {"top": 1, "bottom": -1}.get(self.direction, 0)

    def resolve_genes(self, pool: GenePool) -> list[GeneHit]:
        """Look up gene detail from the shared pool."""
        return [pool[gid] for gid in self.gene_ids if gid in pool]

    def mean_input_value(self, pool: GenePool) -> float:
        """Mean ranking score across gene hits — serves as pseudo-NES.

        Requires the gene pool because gene data is not embedded.
        """
        hits = self.resolve_genes(pool)
        if not hits:
            return 0.0
        return sum(h.input_value for h in hits) / len(hits)
```

### 1.3 `CategoryGSEA` — all terms for one category in one contrast

This is analogous to clusterProfiler's `enrichResult`. It's the natural unit for per-category visualizations: ridge plots, upset plots, term-network graphs. It owns the shared gene pool for its contrast (or receives it from the parent).

```python
class CategoryGSEA(BaseModel):
    """Enrichment results for one gene set category in one contrast.

    Analogous to clusterProfiler's enrichResult — the natural unit
    for single-category visualizations (ridge plot, upset plot, term network).

    Holds a reference to the shared gene pool so that terms can
    resolve their gene detail without embedding copies.
    """

    contrast: str                 # e.g., "Therapy1_vs_Control"
    category: str                 # e.g., "GO Process", "KEGG", "Monarch"
    terms: list[TermGSEA]
    gene_pool: GenePool           # shared across all terms in this category

    @property
    def n_significant(self) -> int:
        return len(self.terms)

    def filter_by_fdr(self, threshold: float = 0.05) -> "CategoryGSEA":
        """Return a new CategoryGSEA with only terms below FDR threshold."""
        return self.model_copy(update={
            "terms": [t for t in self.terms if t.false_discovery_rate < threshold]
        })

    def top_terms(self, n: int = 100) -> "CategoryGSEA":
        """Return top N terms by FDR."""
        sorted_terms = sorted(self.terms, key=lambda t: t.false_discovery_rate)
        return self.model_copy(update={"terms": sorted_terms[:n]})

    def to_polars(self) -> pl.DataFrame:
        """Flatten to a Polars DataFrame (one row per term, no gene explosion).

        Includes computed columns: gene_ratio, mean_input_value, direction_nr.
        """
        ...

    def to_exploded_polars(self) -> pl.DataFrame:
        """Flatten to a Polars DataFrame (one row per gene hit per term).

        Resolves gene_ids through the gene_pool to produce the
        exploded format that some seaborn-based plots need.
        """
        ...
```

**Note on gene pool ownership:** The gene pool is logically per-contrast (a gene has the same input_value in every category). When constructing `CategoryGSEA`, the pool can be either (a) a filtered subset containing only genes in this category's terms, or (b) the full contrast-level pool shared by reference. Option (b) is simpler and uses less memory (Python dicts are shared by reference, not copied). The `GSEAResult` constructor builds one pool per contrast and passes it to all `CategoryGSEA` instances for that contrast.

### 1.4 Container classes — for multi-category and multi-contrast views

```python
class MultiCategoryGSEA(BaseModel):
    """Multiple categories for one contrast.

    Use case: overview visualizations, heatmaps across categories.
    All CategoryGSEA instances share the same gene pool (same contrast).
    """

    contrast: str
    categories: dict[str, CategoryGSEA]     # keyed by category name

    @classmethod
    def from_categories(cls, contrast: str, cats: list[CategoryGSEA]) -> "MultiCategoryGSEA":
        return cls(contrast=contrast, categories={c.category: c for c in cats})

    def __getitem__(self, category: str) -> CategoryGSEA:
        return self.categories[category]

    @property
    def category_names(self) -> list[str]:
        return list(self.categories.keys())


class MultiContrastGSEA(BaseModel):
    """Same category across multiple contrasts.

    Use case: cross-contrast comparison plots (dotplot, upset).
    Each CategoryGSEA has its own gene pool (different contrasts
    have different input values for the same gene).
    """

    category: str
    contrasts: dict[str, CategoryGSEA]      # keyed by contrast name

    @classmethod
    def from_categories(cls, category: str, cats: list[CategoryGSEA]) -> "MultiContrastGSEA":
        return cls(category=category, contrasts={c.contrast: c for c in cats})

    def __getitem__(self, contrast: str) -> CategoryGSEA:
        return self.contrasts[contrast]

    @property
    def contrast_names(self) -> list[str]:
        return list(self.contrasts.keys())

    @property
    def all_term_ids(self) -> set[str]:
        """Union of all term IDs across contrasts."""
        return {t.term_id for c in self.contrasts.values() for t in c.terms}
```

### 1.5 Top-level result container

```python
class GSEAResult(BaseModel):
    """Complete GSEA result for a workunit — all contrasts, all categories.

    This is what you get after loading a long XLSX or collecting API responses.
    Slice it by category or contrast to get the visualization-ready containers.

    Owns the gene pools (one per contrast). Each CategoryGSEA receives
    its contrast's pool by reference during construction.
    """

    workunit_id: str
    species: int
    category_names: list[str]               # available category names
    contrast_names: list[str]               # available contrast names

    # Gene pools: one per contrast (a gene has the same score in every category)
    gene_pools: dict[str, GenePool]         # keyed by contrast name

    # Data: nested dict for natural access — data[contrast][category]
    data: dict[str, dict[str, CategoryGSEA]]

    def get_category(self, contrast: str, category: str) -> CategoryGSEA:
        """Get one category for one contrast."""
        return self.data[contrast][category]

    def get_multi_category(self, contrast: str) -> MultiCategoryGSEA:
        """All categories for one contrast → for overview plots."""
        cats = list(self.data[contrast].values())
        return MultiCategoryGSEA.from_categories(contrast, cats)

    def get_multi_contrast(self, category: str) -> MultiContrastGSEA:
        """One category across all contrasts → for comparison plots."""
        cats = [self.data[c][category] for c in self.contrast_names if category in self.data[c]]
        return MultiContrastGSEA.from_categories(category, cats)

    @classmethod
    def from_long_xlsx(cls, path: Path, workunit_id: str, species: int) -> "GSEAResult":
        """Parse the long XLSX into the full object model.

        Construction flow:
        1. Read XLSX into Polars DataFrame
        2. For each contrast: parse comma-separated protein columns into
           a GenePool (one GeneHit per unique protein_id)
        3. For each (contrast, category): build TermGSEA objects with
           gene_ids referencing the contrast's pool
        4. Wrap in CategoryGSEA with the shared pool
        """
        ...

    @classmethod
    def from_long_polars(cls, df: pl.DataFrame, workunit_id: str, species: int) -> "GSEAResult":
        """Parse a long-format Polars DataFrame into the full object model."""
        ...
```

---

## How This Maps to Visualizations

| Visualization | Current input | Proposed input |
|---------------|---------------|----------------|
| Ridge plot (one category) | `pl.DataFrame` filtered by contrast+category | `CategoryGSEA` |
| Upset plot (one category) | `pl.DataFrame` filtered by contrast+category | `CategoryGSEA` |
| Term-term network | `pl.DataFrame` filtered by contrast+category | `CategoryGSEA` |
| Enrichment dotplot (cross-contrast) | `pl.DataFrame` with multiple contrasts | `MultiContrastGSEA` |
| Upset across contrasts | `pl.DataFrame` with multiple contrasts | `MultiContrastGSEA` |
| Overview heatmap | `pl.DataFrame` grouped by category | `MultiCategoryGSEA` |
| Per-contrast summary table | Full `pl.DataFrame` | `GSEAResult.get_multi_category(contrast)` |

The **key win**: plot functions declare their input type instead of silently expecting specific DataFrame columns. A type error at the call site is much easier to debug than a `ColumnNotFoundError` deep inside a plotting function.

---

## Phase 2: ORA Equivalents

ORA has no per-gene scores, so the model is simpler:

```python
class ORATermResult(BaseModel):
    """A single term from STRING-DB ORA."""

    term_id: str
    description: str
    category: str
    number_of_genes: int              # overlap count
    number_of_genes_in_background: int | None = None
    p_value: float
    fdr: float
    input_genes: list[str] = []       # parsed from comma-separated

    @property
    def gene_ratio(self) -> float:
        if not self.number_of_genes_in_background:
            return 0.0
        return self.number_of_genes / self.number_of_genes_in_background


class CategoryORA(BaseModel):
    """ORA results for one gene set category."""

    category: str
    terms: list[ORATermResult]
```

ORA containers (`MultiCategoryORA`, etc.) can follow the same pattern if needed.

---

## Phase 3: Serializers

Methods on the model classes + standalone functions for interchange formats:

```python
# On GSEAResult:
def to_long_xlsx(self, path: Path) -> None: ...
def to_pivoted_xlsx(self, path: Path) -> None: ...
def to_merged_xlsx(self, path: Path) -> None: ...

# Standalone for interoperability:
def to_clusterprofiler_df(cat: CategoryGSEA) -> pl.DataFrame:
    """Export as clusterProfiler-compatible DataFrame."""
    ...

def to_gseapy_df(cat: CategoryGSEA) -> pl.DataFrame:
    """Export as GSEApy-compatible DataFrame."""
    ...
```

### Replacing `gsea_result_processor.py`

The current `_results_to_dataframe()` → pivot → merge flow would become:

```python
result = GSEAResult.from_long_xlsx(path, workunit_id, species)
result.to_long_xlsx(out / "long.xlsx")
result.to_pivoted_xlsx(out / "pivoted.xlsx")
result.to_merged_xlsx(out / "merged.xlsx")
```

---

## Phase 4: Column Name Mapping Reference

For documentation and potential export. Maps STRING-DB names to community equivalents:

| STRING-DB GSEA | clusterProfiler | GSEApy | fgsea | g:Profiler |
|----------------|-----------------|--------|-------|------------|
| `termID` | `ID` | `Term` | `pathway` | `term_id` |
| `termDescription` | `Description` | — | — | `term_name` |
| `category` | — | `Gene_set` | — | `source` |
| `enrichmentScore` | `enrichmentScore` | `ES` | `ES` | — |
| `meanInputValue` (computed) | `NES` (approximate) | `NES` | `NES` | — |
| `falseDiscoveryRate` | `p.adjust` | `FDR q-val` | `padj` | `p_value` |
| — (no raw p) | `pvalue` | `NOM p-val` | `pval` | — |
| `genesInSet` | `setSize` | — | `size` | `term_size` |
| `genesMapped` | `Count` | — | — | `overlap_size` |
| `direction` | sign of NES | sign of NES | sign of NES | — |
| `proteinLabels` | `core_enrichment` | `Lead_genes` | `leadingEdge` | `intersection` |

---

## Data Flow: Current vs Proposed

### Current

```
STRING API TSV → pd.read_csv → concat all contrasts → add directionNR/num_contrasts
    → write long XLSX → pivot → write pivoted/merged XLSX
    → read long XLSX → filter_by_FDR → add_gene_ratio → explode_protein_columns
    → summarize_terms → pass DataFrame to 15 × 5 plot functions
```

Every stage passes `pl.DataFrame`. Columns are implicit. The `explode_protein_columns` step blows up row count ~56x.

### Proposed

```
STRING API TSV
    → parse into GSEAResult:
        1. Build GenePool per contrast (one GeneHit per unique protein_id)
        2. Build TermGSEA with gene_ids referencing the pool
        3. Group into CategoryGSEA (each holds ref to contrast's pool)
    → GSEAResult.to_long_xlsx() / to_pivoted_xlsx() / to_merged_xlsx()
    → GSEAResult.from_long_xlsx() (for report generation)
    → result.get_category("contrast", "GO Process") → CategoryGSEA → plot_ridges(cat)
    → result.get_multi_contrast("GO Process") → MultiContrastGSEA → plot_dotplot(mc)
    → result.get_multi_category("contrast") → MultiCategoryGSEA → plot_heatmap(mc)
```

Gene data lives in a shared pool per contrast — no duplication, no explosion. Terms reference genes by ID. When a plot function does need a flat DataFrame (e.g., for seaborn), `CategoryGSEA.to_polars()` produces term-level rows, and `.to_exploded_polars()` resolves gene IDs through the pool to produce the exploded format.

---

## Priority and Sequencing

| Phase | What | Depends on | Effort |
|-------|------|------------|--------|
| 1 | Pydantic models (GeneHit → TermGSEA → CategoryGSEA → containers → GSEAResult) | Nothing | Medium |
| 2 | ORA equivalents (ORATermResult → CategoryORA) | Nothing | Small |
| 3 | Serializers (XLSX writers, clusterProfiler export) | Phase 1 | Medium |
| 4 | Parser: `GSEAResult.from_long_xlsx()` + `from_long_polars()` | Phase 1 | Medium |
| 5 | Migrate plotting functions to accept typed containers | Phase 1 + TODO_improve_plotting_api.md | Large |

Phase 1 can be built and tested standalone. Phase 4 (the parser) is the integration point — it replaces `_results_to_dataframe()`. Phase 5 is the big payoff but overlaps with the plotting API cleanup (TODO_improve_plotting_api.md Phase D).

---

## Open Questions

1. **pydantic vs dataclasses?** Pydantic gives validation, aliases, JSON serialization, and `model_copy()` for free. Cost: one new dependency. Recommendation: **pydantic** — this project is not dependency-shy (already has cyclopts, polars, matplotlib, networkx, etc.).

2. **Gene pool granularity.** The pool is logically per-contrast (a gene's input_value and rank are the same regardless of which category matched it). All `CategoryGSEA` instances for the same contrast share the same pool by reference. Cross-contrast containers (`MultiContrastGSEA`) hold categories with different pools — this is correct since the same gene can have different scores in different contrasts.

3. **Comma-separated protein column parsing.** The `§COMMA§` protection logic (for gene symbols containing commas followed by 1-2 digits) needs to be preserved in the `GSEAResult` constructor when building the gene pool. This is a one-time parse per contrast, replacing the current `explode_protein_columns()` which re-parses per operation.

4. **Should `CategoryGSEA` hold filtered or unfiltered terms?** Options: (a) always store all terms, filter via methods, (b) filter at construction time. Recommendation: **(a)** — store all, provide `filter_by_fdr()` and `top_terms()` methods that return new instances. This preserves the raw data for summary statistics ("X of Y terms significant").

5. **Should the `meanInputValue` pseudo-NES replace `directionNR`?** The mean input value is strictly more informative (continuous, signed). `directionNR` (1/-1/0) is just a discretized version. We could keep both: `direction_nr` for backward compatibility, `mean_input_value` as the richer alternative. Plotting code can choose which to use.

6. **Gene pool serialization.** When writing to XLSX, the gene pool needs to be flattened back into comma-separated protein columns (for backward compatibility with existing consumers). When serializing the object model itself (e.g., to JSON for caching), the pool is serialized once and terms store only IDs — much more compact than the current long XLSX with duplicated protein strings.
