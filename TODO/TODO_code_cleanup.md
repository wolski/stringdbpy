# Code Cleanup TODO — stringdbpy

Date: 2026-03-30

## Context

This review traces all reachable code from the formal entry points:

- **CLI entry points** (pyproject.toml): `string_gsea_run`, `string_gsea_render_reports`, `string_gsea_render_marimo`, `string_ora_run`, `string_gsea_write_config`
- **Snakefile** (`tests/Snakefile`): calls `string_gsea_run` + `string_gsea_render_marimo`
- **External consumer A383** (`slurmworker/config/A383_STRING_GSEA_HTML/Snakefile`): calls `string_gsea_run` + `string_gsea_render_reports`
- **External consumer A373** (`slurmworker/config/A373_STRING_GSEA/run.sh`): calls `string_gsea_bfabric` (see Phase 2)

---

## Phase 0 — Stabilize Tests First

Before any code deletion, the test suite must be reliable so we can validate each step.

### 0.1 Current Test Classification

| Test File | Type | Notes |
|-----------|------|-------|
| `test_taxon_utils.py` | fast/unit | Conditional `@skipif` on local data files |
| `test_gsea_utilities.py` | fast/unit | Local file operations only |
| `test_gsea_session.py` | fast/unit | YAML serialization round-trips |
| `test_render_reports.py` | fast/unit | Fully mocked subprocess |
| `test_prep_data_input.py` | fast/unit | tmp_path fixtures |
| `test_network.py` | fast/unit | Local DataFrames, `pytest.skip` on missing deps |
| `test_diff_xlsx_integration.py` | fast/unit | Reads test XLSX (name is misleading) |
| **`test_get_species.py`** | **integration** | **LIVE STRING API calls, unconditionally** |
| `test_config_integration.py` | fast/unit | Fully mocked despite name |
| `test_gsea_results.py` | integration | Guarded by `RUN_STRING_GSEA_INTEGRATION=1` env var |

**Problems:**
1. No pytest markers (`@pytest.mark.slow`, `@pytest.mark.integration`) to separate fast from slow.
2. `test_get_species.py` hits the STRING API unconditionally — will fail offline or slow down CI.
3. Test file names are misleading (`test_diff_xlsx_integration.py` is actually unit, `test_config_integration.py` is mocked).

### 0.2 Action: Add pytest markers

- Define markers in `pyproject.toml` or `conftest.py`: `integration` (hits STRING API), `ci` (fast, <30s).
- Mark `test_get_species.py` as `@pytest.mark.integration`.
- Default `pytest` invocation should skip `integration` tests.
- Snakefile `ci` rule already tags small datasets (`human_rnk_2848501`, `mouse_xlsx`, `yeast_rnk` — ~20MB total).

### 0.3 Action: Finalize test data layout

**Current state (incomplete migration):**
- Old files deleted from working tree but not committed:
  ```
  D tests/data/2848501.zip
  D tests/data/DE_mouse_fasta_xlsx.zip
  D tests/data/DE_yeast_fasta_rnk.zip
  D tests/data/dummy_d/*
  D tests/data/dummy_out/*
  D tests/data/no_matching_zip/*
  ```
- New layout: `tests/data/datasets/` (9 datasets with `params.yml`) and `tests/data/fixtures/` (3 fixture dirs).

**Test datasets inventory:**

| Dataset | Species | Input | Size | Tags | CI? |
|---------|---------|-------|------|------|-----|
| `human_rnk_2848501` | Human | RNK | 8.4 MB | small | **yes** |
| `mouse_xlsx` | Mouse | XLSX | 9.5 MB | small | **yes** |
| `yeast_rnk` | Yeast | RNK | 2.0 MB | small | **yes** |
| `human_rnk_2881688` | Human | RNK | 76 MB | large | no |
| `human_rnk_2958180_single` | Human | RNK | 211 MB | large, single | no |
| `human_xlsx_20250813` | Human | XLSX | 44 MB | medium | no |
| `human_xlsx_20250908` | Human | XLSX | 39 MB | medium | no |
| `human_xlsx_primary_58167` | Human | XLSX | 80 MB | large | no |
| `human_xlsx_secondary_58167` | Human | XLSX | 82 MB | large | no |

**Actions:**
- Commit the old test data deletions + new layout.
- Add `tests/.snakemake/` to `.gitignore`.
- Decide on `tests/errors/` (~355 MB of debugging artifacts: `2990472.zip`, `GSEA_STRING_GSEA_HTML_334883_2989749/`, `o33436_Hypoxia_vs_Normoxia_IMPUTE/`, `testing/`). These contain real rnk/xlsx data but are unstructured. Either delete or promote useful ones to `tests/data/datasets/` with proper `params.yml`.

### 0.4 Action: Rename misleading test files

- `test_diff_xlsx_integration.py` → `test_diff_xlsx.py` (it's a unit test)
- `test_config_integration.py` → `test_config.py` (it's fully mocked)

---

## Phase 0.5 — Configurable API Base URL (GitHub issue #13)

The STRING-DB API base URL `https://version-12-0.string-db.org/api` is hardcoded in 5 locations.
Once tests are stable (Phase 0), this is the first code change — it centralizes the URL
and unblocks integration test configuration.

### Hardcoded locations

| File | Location | Usage |
|------|----------|-------|
| `gsea_config.py:84` | `_fetch_api_key()` default param | Config setup |
| `gsea_session.py:22-23` | Class attr `end_point_status` | Polling URL |
| `string_gsea_builder.py:48` | Local in `_submit_single()` | Job submission |
| `string_ora_run.py:15` | Module constant `STRING_API_BASE` | ORA endpoints |
| `get_species.py:34` | Local in `_fetch_ncbi_taxon_ids()` | Species lookup |

### Actions

1. Add `api_base_url: str = "https://version-12-0.string-db.org/api"` to `GSEAConfig` (optional with default — existing config.toml files keep working).
2. Update `read_toml()` and `from_dict()` to handle the new field.
3. Replace all 5 hardcoded URLs with references to `config.api_base_url`.
4. `get_species_taxon()` gains an `api_base_url` parameter (call site in `string_gsea_run.py` has config).
5. Session YAML includes the URL automatically via `asdict()`.

---

## Phase 1 — Define Entry Points

### 1.1 `string_gsea_bfabric` — removed but still needed

**Finding:** `string_gsea_bfabric` was removed from stringdbpy in commit `ba7c0b5` (Sept 2025). However:
- **A373** (`slurmworker/config/A373_STRING_GSEA/`) has its own `pyproject.toml` that **still defines the entry point locally**: `"string_gsea_bfabric" = "string_gsea.scripts.string_gsea_bfabric:app"`. So A373 works from its locked/vendored copy, but only as long as that lock file isn't regenerated.
- **A383** (`slurmworker/config/A383_STRING_GSEA_HTML/`): `run.sh` calls `string_gsea_bfabric` but is **dead code** — `app.yml` invokes snakemake directly (`command: snakemake --cores 1 -s .../Snakefile -d`), never referencing `run.sh`. The Snakefile calls `string_gsea_run` and `string_gsea_render_reports`. `run.sh` can be deleted from A383.

**Decision needed:** A373 runs GSEA without HTML reporting. Options:
1. Restore `string_gsea_bfabric` as an entry point (or alias to `string_gsea_run`).
2. Create a Snakefile for A373 (like A383 has) and update `run.sh` to use `string_gsea_run` directly.
3. Unify A373 and A383 into a single Snakefile-based config with a flag for report generation.

### 1.2 Current entry points (keep)

```
string_gsea_run          → scripts/string_gsea_run.py:app       (main GSEA workflow)
string_gsea_render_reports → scripts/render_reports.py:app       (Quarto reports)
string_gsea_render_marimo  → scripts/render_marimo_reports.py:app (Marimo reports)
string_ora_run           → scripts/string_ora_run.py:app         (ORA analysis)
string_gsea_write_config → scripts/write_config.py:app           (config setup)
```

---

## Phase 2 — Safe Deletes (no behavior change)

### 2.1 Merge `config.py` into `gsea_config.py`

**Finding:** The two files are nearly identical, but `gsea_config.py` is the **more advanced version**:
- Has `required` as a **class attribute** (reusable) vs. a local variable in `config.py`
- Has an extra `from_dict()` classmethod for constructing config from dicts
- Is the version imported by all production code

`config.py` is the older version. `gsea_config.py` imports from it on line 11 (`from string_gsea.config import GSEAConfig`) but immediately shadows that import with its own `GSEAConfig` redefinition.

**Action:** Delete `config.py`, remove the dead import on line 11 of `gsea_config.py`. No behavior change — `gsea_config.py` already has all functionality.

### 2.2 Remove `find_zip_files()` from `gsea_utilities.py`

- **Location**: `src/string_gsea/gsea_utilities.py`, around line 28
- **Why**: Never called from any entry point. Only referenced by `tests/test_gsea_utilities.py`.
- **Action**: Remove the function and its corresponding test.

### 2.3 Remove `test_run()` from `string_gsea_run.py`

- **Location**: `src/string_gsea/scripts/string_gsea_run.py`, around line 93
- **Why**: Debug function. Commented out in `__main__` block. Never called.
- **Action**: Delete the function.

### 2.4 Remove `get_session_file_path()` from `_report_utils.py`

- **Location**: `src/string_gsea/marimo_reports/_report_utils.py`, line 82
- **Why**: Defined but never called. **Checked Quarto integration:** no references in any `.qmd` file either. Completely orphaned.
- **Action**: Delete the function.

### 2.5 Delete `src/string_gsea/docs/index.qmd.bak`

- **Why**: Orphaned backup file, untracked by git. Not referenced by any `_quarto.yml` profile.
- **Action**: Delete.

---

## Phase 3 — After confirming CircularGraph.qmd is not needed

### 3.1 Delete `src/string_gsea/docs/python_notebooks/CircularGraph.qmd`

- **Why**: Not referenced by any `_quarto.yml` profile or render script. Development/exploration notebook outside production pipeline.
- **Action**: Delete the file and the `python_notebooks/` directory if empty.

### 3.2 Remove unused visualization functions from `network.py`

These functions are **only** reachable from the orphaned `CircularGraph.qmd` and from tests. Production Quarto templates and marimo reports use a different subset of `network.py` (`filter_by_FDR`, `add_gene_ratio`, `explode_protein_columns`, `summarize_terms`, `TermNetworkBuilder`, `TermNetworkPlotter`).

Functions to remove from `src/string_gsea/network.py`:

| Function | Notes |
|----------|-------|
| `bipartite_barycenter_layout()` | Never called anywhere. Superseded by `bipartite_hybrid_layout()`. |
| `bipartite_hybrid_layout()` | Only called by `plot_network_graph_plotly()`, which is itself dead. |
| `make_network()` | Only called by `make_network_with_colors()`, which is dead. |
| `make_network_with_colors()` | Only used by `CircularGraph.qmd` and tests. |
| `plot_network_graph()` | Only used by `CircularGraph.qmd` and tests. |
| `interactive_cytoscape()` | Only used by `CircularGraph.qmd` and tests. |
| `plot_network_graph_plotly()` | Only used by `CircularGraph.qmd` and tests. |
| `build_tooltip()` | Only called by `interactive_cytoscape()` and `plot_network_graph_plotly()`. |

- **Action**: Remove all of the above. Update `test_network.py` to only test functions used by production code paths.

---

## Phase 4 — Minor Cleanup

### 4.1 `ReportConfig.from_env()` in `_report_utils.py`

- Keep but add a docstring clarifying it's a dev/interactive convenience path. Production uses `from_args()`.

### 4.2 `tests/.snakemake/` directory

- Add to `.gitignore`.

---

## Summary

| Phase | Items | Risk | Prerequisite |
|-------|-------|------|-------------|
| Phase 0 | Test stabilization (markers, data layout, renames) | None | — |
| Phase 0.5 | Configurable API base URL (issue #13) | None | Phase 0 done |
| Phase 1 | Entry point decisions (A373 `string_gsea_bfabric`) | Medium | Needs decision |
| Phase 2 | 5 safe deletes (dead code with no callers) | None | Phase 0 done |
| Phase 3 | CircularGraph.qmd + ~8 dead functions in network.py | Low | Phase 0 done |
| Phase 4 | Minor housekeeping | None | Any time |

**Biggest concern**: `network.py` (~600 lines) where roughly half the public API is unreachable from production code paths. This is the main legacy of the stopped refactoring.
