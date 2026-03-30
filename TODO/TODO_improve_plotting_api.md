# TODO: Improve Plotting API

Date: 2026-03-30

## Root Problem

The mouse_xlsx CI dataset (single contrast) generates a **55MB** marimo HTML report that hangs the export. Root cause: **no limit on categories or terms** rendered. The "Publications" category alone has 7,472 rows after FDR filtering. The report blindly iterates all 15 categories x 5 plot types = 75 plots, most of them enormous.

Additionally, the plotting code is scattered across 6 files with no consistent API, no input validation, and several global state mutations that cause side effects.

---

## Current Plotting Files

| File | Functions | Purpose | Used By |
|------|-----------|---------|---------|
| `gsea_plotting.py` | `plot_single_ridge`, `plot_term_ridges`, `make_upset`, `make_upset_contrasts_terms` | Ridge plots, upset plots | marimo reports, QMD |
| `dotplot_endrichment.py` | `prepare_data_for_plotting`, `plot_enrichment_scatter`, `add_custom_legends`, `dotplot_enrichment` | Enrichment dotplot | marimo report_multiple, QMD |
| `network.py` | `filter_by_FDR`, `add_gene_ratio`, `explode_protein_columns`, `summarize_terms`, `make_network`, `assign_node_sizes`, `assign_node_colors`, `make_network_with_colors`, `plot_network_graph`, `interactive_cytoscape`, `plot_network_graph_plotly`, `build_tooltip`, `bipartite_hybrid_layout`, `bipartite_barycenter_layout` | Data transforms + bipartite network viz | marimo reports, QMD, tests |
| `TermNetworkBuilder.py` | `build_shared_counts`, `build_contrast_counts`, `compute_node_sizes` | Term-term network data | `TermNetworkPlotter` |
| `TermNetworkPlotter.py` | `plot_network`, `draw_panel`, `compute_full_layout`, `draw_legend_panel`, `get_figure_legend` | Term-term network viz | marimo report_networks, QMD |
| `cluster_genesets.py` | `pivot_to_wide`, `make_nested_dict`, `convert_to_binary`, `plot_term_distance_heatmap` | Heatmap + clustering | marimo report_networks, QMD |

### Problems

1. **No term/category limits** — renders all categories including Publications (7,472 rows after FDR filter)
2. **Data transforms mixed with plotting** — `filter_by_FDR`, `add_gene_ratio`, `explode_protein_columns`, `summarize_terms` are in `network.py` but have nothing to do with network visualization
3. **Global state mutations**:
   - `TermNetworkPlotter.plot_network()` sets `plt.rcParams["figure.dpi"] = 300` (never reset, affects all subsequent figures)
   - `plot_term_ridges()` calls `sns.despine(left=True, bottom=True)` globally
4. **Inconsistent return types** — some functions return Figure, some return None + call `plt.show()`, some return UpSet objects
5. **No input dataclass** — each function expects specific columns but there's no schema definition
6. **Dead code in network.py** — `make_network`, `make_network_with_colors`, `plot_network_graph`, `interactive_cytoscape`, `plot_network_graph_plotly`, `build_tooltip`, `bipartite_hybrid_layout`, `bipartite_barycenter_layout` (see TODO_code_cleanup.md Phase 3)
7. **Typo in filename** — `dotplot_endrichment.py` (should be `dotplot_enrichment.py`)

---

## Proposed Structure

```
src/string_gsea/
├── data_transforms.py          # filter_by_FDR, add_gene_ratio, explode_protein_columns, summarize_terms
├── plots/
│   ├── __init__.py             # Public API re-exports
│   ├── _types.py               # Input dataclass(es)
│   ├── ridgeplot.py            # plot_term_ridges (+ plot_single_ridge as private)
│   ├── upset.py                # make_upset, make_upset_contrasts_terms
│   ├── dotplot.py              # dotplot_enrichment (+ helpers as private)
│   ├── heatmap.py              # plot_term_distance_heatmap (from cluster_genesets.py)
│   └── term_network.py         # TermNetworkBuilder + TermNetworkPlotter (merged)
├── cluster_genesets.py         # pivot_to_wide, make_nested_dict, convert_to_binary (data utils, no plotting)
└── network.py                  # DELETED after moving transforms + removing dead code
```

---

## Phase A — Immediate Fix: Add top-N term filter (unblocks CI)

The report hangs because it tries to render thousands of terms per category. Add a `max_terms` parameter to the data pipeline that keeps only the top N terms by FDR per category.

### A.1 Add `select_top_terms()` function

Add to the data transform pipeline (currently in `network.py`, later in `data_transforms.py`):

```python
def select_top_terms(df: pl.DataFrame, max_terms: int = 100) -> pl.DataFrame:
    """Keep only the top max_terms terms per (contrast, category) by FDR."""
    return (
        df.sort(["contrast", "category", "falseDiscoveryRate"])
        .with_columns(
            pl.col("termID")
            .rank(method="dense")
            .over(["contrast", "category"])
            .alias("_rank")
        )
        .filter(pl.col("_rank") <= max_terms)
        .drop("_rank")
    )
```

### A.2 Integrate into marimo reports

In `report_networks.py` and `report_multiple.py`, call after `filter_by_FDR`:

```python
df = filter_by_FDR(df, config.fdr_threshold, config.genes_mapped_threshold)
df = select_top_terms(df, max_terms=100)  # NEW
```

### A.3 Add `--max-terms` CLI parameter

In `render_marimo_reports.py` and `render_reports.py`, add parameter (default 100).

---

## Phase B — Fix global state mutations

### B.1 Remove `plt.rcParams["figure.dpi"] = 300` from `TermNetworkPlotter.plot_network()`

**File:** `TermNetworkPlotter.py`, line 259

This sets DPI globally and is never reset. The report already saves PNGs at `dpi=100` in `plt.savefig()`. Remove the global mutation. If high-res is needed, pass `dpi` to `savefig()` at the call site.

### B.2 Remove global `sns.despine()` from `plot_term_ridges()`

**File:** `gsea_plotting.py`

`sns.despine(left=True, bottom=True)` affects all subsequent axes. Move despine inside the per-axes loop or scope it to the figure.

### B.3 Consistent return types

All plot functions should **return the figure** and never call `plt.show()`. The caller (marimo report, QMD, interactive session) decides when to show.

Functions that currently call `plt.show()` instead of returning:
- `plot_network()` in `TermNetworkPlotter.py` — calls `plt.show()` at end
- `dotplot_enrichment()` in `dotplot_endrichment.py` — calls `plt.show()`
- `make_upset_contrasts_terms()` in `gsea_plotting.py` — calls `plt.show()`

---

## Phase C — Define input dataclass

### C.1 Standard enrichment DataFrame schema

The pipeline expects a Polars DataFrame with specific columns at each stage. Define this as a dataclass or TypedDict for documentation and validation:

```python
@dataclass
class EnrichmentData:
    """Standard GSEA enrichment result schema.

    Required columns:
        contrast, category, termID, termDescription,
        genesMapped, genesInSet, enrichmentScore, direction,
        falseDiscoveryRate, proteinIDs, proteinLabels,
        proteinInputLabels, proteinInputValues, proteinRanks

    Added by pipeline:
        directionNR (int): -1/0/1 for bottom/unknown/top
        geneRatio (float): genesMapped / genesInSet
        meanInputValues (float): mean protein value per term
    """
    df: pl.DataFrame

    def validate(self) -> None:
        """Check required columns exist."""
        ...
```

This is mostly for documentation — Polars doesn't enforce schemas at type-check time.

---

## Phase D — Reorganize into `plots/` module

### D.1 Extract data transforms from `network.py`

Move to `data_transforms.py`:
- `filter_by_FDR()`
- `add_gene_ratio()`
- `explode_protein_columns()`
- `summarize_terms()`
- `select_top_terms()` (new from Phase A)

### D.2 Create `plots/` subpackage

Move plotting functions as described in "Proposed Structure" above. Keep backward-compatible re-exports in `__init__.py` temporarily if needed by QMD templates (they import from `string_gsea.network`).

### D.3 Merge `TermNetworkBuilder` + `TermNetworkPlotter`

These are always used together and the builder has no other consumers. Merge into `plots/term_network.py`.

### D.4 Delete dead code from `network.py`

After extracting transforms and live plotting code, delete:
- `make_network()`, `assign_node_sizes()`, `assign_node_colors()`, `make_network_with_colors()`
- `plot_network_graph()`, `interactive_cytoscape()`, `plot_network_graph_plotly()`
- `build_tooltip()`, `bipartite_hybrid_layout()`, `bipartite_barycenter_layout()`

(Already tracked in `TODO_code_cleanup.md` Phase 3)

### D.5 Fix filename typo

`dotplot_endrichment.py` → `plots/dotplot.py`

### D.6 Move `pivot_to_wide`, `make_nested_dict`, `convert_to_binary` out of `cluster_genesets.py`

These are data utilities, not clustering. Move to `data_transforms.py`. Keep `plot_term_distance_heatmap` → `plots/heatmap.py`.

---

## Data Flow (current vs proposed)

### Current

```
XLSX → filter_by_FDR (network.py)
     → add_gene_ratio (network.py)
     → explode_protein_columns (network.py)
     → summarize_terms (network.py)
     → make_nested_dict (cluster_genesets.py)
     → 15 categories × 5 plots each = up to 75 renders (NO LIMIT)
```

### Proposed

```
XLSX → filter_by_FDR (data_transforms.py)
     → select_top_terms(max_terms=100) (data_transforms.py)    ← NEW
     → add_gene_ratio (data_transforms.py)
     → explode_protein_columns (data_transforms.py)
     → summarize_terms (data_transforms.py)
     → make_nested_dict (data_transforms.py)
     → ≤15 categories × ≤100 terms × 5 plots = bounded renders
```

---

## Priority

| Phase | What | Urgency | Effort |
|-------|------|---------|--------|
| A | Top-N term filter | **High** — unblocks CI | Small |
| B | Fix global state mutations | **High** — causes side effects | Small |
| C | Input dataclass | Low — documentation | Small |
| D | Reorganize into plots/ | Low — cleanup | Medium |
