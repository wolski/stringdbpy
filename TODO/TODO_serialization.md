# TODO: GSEAResult Serialization & Export

Date: 2026-03-31

## Motivation

`GSEAResult` is now the typed source of truth in the pipeline, but we need ways to persist, exchange, and consume it from other languages (especially R).

---

## Phase 1: XLSX Export — DONE (via `to_polars_long()`)

`GSEAResult.to_polars_long()` returns a Polars DataFrame. The caller writes XLSX however they want. `GSEAResultProcessor` now uses this internally. No XLSX writer on the model itself — that's the caller's concern.

---

## Phase 2: JSON Serialization — DONE

`to_json(path)` / `from_json(path)` on `GSEAResult`. Also `to_dict()` / `from_dict()` on every class in the hierarchy.

Design decisions:
- Gene pool stored once per contrast in JSON (inside `MultiCategoryGSEA`), not per category.
- On deserialization, one `GenePool` rebuilt per contrast and shared by reference across all categories.
- Manual `to_dict()` / `from_dict()` (~60 lines) instead of pydantic — avoids new dependency, keeps shared pool pattern clean.
- `tuple` fields (`gene_ids`, `terms`) serialize as JSON arrays, restored as tuples on load.

### JSON structure

```json
{
  "data": {
    "<contrast_name>": {
      "contrast": "...",
      "gene_pool": { "<protein_id>": {"protein_id": "...", "label": "...", ...}, ... },
      "categories": {
        "<category_name>": {
          "category": "...",
          "contrast": "...",
          "terms": [ {"term_id": "...", "gene_ids": [...], ...}, ... ]
        }
      }
    }
  },
  "rank_lists": {
    "<contrast_name>": { "contrast": "...", "entries": {"<gene>": 1.23, ...} }
  }
}
```

---

## Phase 3: R Package `stringGSEAplot` — DONE

R package at `stringGSEAplot/` reads JSON (from Phase 2) and constructs clusterProfiler's `enrichResult` S4 objects for use with enrichplot.

### Package contents

- `read_gsea_json(json_path, categories = NULL)` — entry point, returns named list of `enrichResult` objects keyed by `"contrast::category"`
- `build_enrichResult(category_data, gene_pool, rank_list)` — constructs one `enrichResult` per contrast × category
- `resolve_gene_ids(gene_ids, gene_pool)` — maps STRING protein IDs to labels

### What this unlocks

enrichplot's full visualization suite with zero custom Python plotting:
`dotplot()`, `cnetplot()`, `ridgeplot()`, `upsetplot()`, `treeplot()`, `emapplot()`

### Implication for Python plotting

The custom Python plotting code (`gsea_plotting.py`, `dotplot_endrichment.py`, `TermNetworkPlotter.py`, etc.) can be deprecated. This intersects with `TODO_improve_plotting_api.md`.

---

## Priority

| Phase | Status | Effort |
|-------|--------|--------|
| 1 — XLSX export | DONE | — |
| 2 — JSON round-trip | DONE | — |
| 3 — R package (JSON → enrichResult) | DONE | — |
