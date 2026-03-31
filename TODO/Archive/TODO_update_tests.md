# Implementation Plan: Phase 0 + Phase 0.5

## Context

The test suite needs stabilization before any cleanup or feature work. Currently:
- `test_get_species.py` makes **live STRING API calls unconditionally** (will fail offline)
- `test_gsea_results.py` has `INTEGRATION_FLAG = True` hardcoded, bypassing the env var guard
- No pytest markers to separate fast from integration tests
- Test file names are misleading (e.g., `test_diff_xlsx_integration.py` is actually a unit test)
- Old test data layout is half-migrated (deleted files not committed)
- `tests/.snakemake/` not in `.gitignore`

After stabilizing tests, we centralize the hardcoded STRING API base URL (issue #13) into `GSEAConfig` — a prerequisite for configurable integration testing.

---

## Phase 0: Stabilize Tests

### Step 0.1 — Add pytest markers to `pyproject.toml`

**File:** `pyproject.toml`

Add `[tool.pytest.ini_options]` section:

```toml
[tool.pytest.ini_options]
markers = [
    "integration: tests that make live STRING-DB API calls (deselect with '-m \"not integration\"')",
]
addopts = "-m 'not integration'"
```

This makes `pytest` skip integration tests by default. Run them explicitly with `pytest -m integration` or `pytest -m ""` for all.

### Step 0.2 — Mark integration tests

**File:** `tests/test_get_species.py`
- Add `@pytest.mark.integration` to `test_GetTaxonID_determine_species` (line 37) — the only test that hits the STRING API.
- `test_get_species_from_oxes` and `test_get_ox_fields` are local/unit — leave unmarked.

**File:** `tests/test_gsea_results.py`
- Remove `INTEGRATION_FLAG = True` override (line 11)
- Remove the `@pytest.mark.skipif` decorator
- Replace with `@pytest.mark.integration`
- This is cleaner: one marker system instead of env vars + skipif.

### Step 0.3 — Rename misleading test files

- `tests/test_diff_xlsx_integration.py` → `tests/test_diff_xlsx.py`
- `tests/test_config_integration.py` → `tests/test_config.py`

Use `git mv` for proper tracking.

### Step 0.4 — Add `tests/.snakemake/` to `.gitignore`

**File:** `.gitignore`

Add line: `tests/.snakemake/`

### Step 0.5 — Commit test data migration

Stage the deleted old test data files and the new layout. This is a git housekeeping commit — no code changes. Files to stage:

```
D tests/data/2848501.zip
D tests/data/DE_mouse_fasta_xlsx.zip
D tests/data/DE_yeast_fasta_rnk.zip
D tests/data/dummy_d/
D tests/data/dummy_out/
D tests/data/no_matching_zip/
```

Plus the new `tests/data/datasets/` and `tests/data/fixtures/` directories (already untracked/added).

### Step 0.6 — Verify

```bash
uv run pytest                          # runs only fast/unit tests
uv run pytest -m integration           # runs only integration tests
uv run pytest -m ""                    # runs all tests
```

All unit tests should pass without network access.

---

## Phase 0.5: Configurable API Base URL (issue #13)

### Step 0.5.1 — Add `api_base_url` field to `GSEAConfig`

**File:** `src/string_gsea/gsea_config.py`

```python
@dataclass
class GSEAConfig:
    api_key: str
    fdr: float
    ge_enrichment_rank_direction: int
    caller_identity: str
    creation_date: Optional[str] = None
    api_base_url: str = "https://version-12-0.string-db.org/api"

    required = ["api_key", "fdr", "ge_enrichment_rank_direction", "caller_identity"]
```

Update `read_toml()` (add after line 42):
```python
api_base_url=data.get("api_base_url", "https://version-12-0.string-db.org/api"),
```

Update `from_dict()` (add after line 60):
```python
api_base_url=data.get("api_base_url", "https://version-12-0.string-db.org/api"),
```

Update `_fetch_api_key()` default parameter (line 84): derive from constant:
```python
DEFAULT_API_BASE = "https://version-12-0.string-db.org/api"

def _fetch_api_key(
    url: str = f"{DEFAULT_API_BASE}/json/get_api_key",
) -> Tuple[str, str]:
```

Actually simpler: just reference `GSEAConfig.api_base_url` default isn't accessible as a class attr on a dataclass. Use a module-level constant instead:

```python
STRING_API_BASE_DEFAULT = "https://version-12-0.string-db.org/api"
```

Then use it in the dataclass default and in `_fetch_api_key`.

### Step 0.5.2 — Update `GSEASession`

**File:** `src/string_gsea/gsea_session.py`

Remove the class variable:
```python
end_point_status = "https://version-12-0.string-db.org/api/json/valuesranks_enrichment_status"
```

Add a property:
```python
@property
def end_point_status(self) -> str:
    return f"{self.config_dict.api_base_url}/json/valuesranks_enrichment_status"
```

Note: `end_point_status` is accessed in `string_gsea_builder.py:81` as `self.session.end_point_status` — the property is a drop-in replacement.

YAML serialization: `api_base_url` will be included automatically via `asdict(self.config_dict)` in `to_yaml()`. Deserialization via `GSEAConfig.from_dict()` will pick it up with the default fallback.

### Step 0.5.3 — Update `StringGSEABuilder`

**File:** `src/string_gsea/string_gsea_builder.py`

In `_submit_single()` (line 48), replace:
```python
base = "https://version-12-0.string-db.org/api"
```
with:
```python
base = self.session.config_dict.api_base_url
```

### Step 0.5.4 — Update `string_ora_run.py`

**File:** `src/string_gsea/scripts/string_ora_run.py`

Remove module constant (line 15):
```python
STRING_API_BASE = "https://version-12-0.string-db.org/api"
```

In `string_ora_run()` function, after `config = get_configuration()`, use `config.api_base_url` instead. Thread it to all functions that need it:

- `map_to_string_ids()` — add `api_base_url` parameter
- `get_string_link()` — add `api_base_url` parameter
- `run_enrichment()` — add `api_base_url` parameter

### Step 0.5.5 — Update `get_species.py`

**File:** `src/string_gsea/get_species.py`

In `GetTaxonID._fetch_ncbi_taxon_ids()` (line 34), add `api_base_url` parameter:
```python
@staticmethod
def _fetch_ncbi_taxon_ids(identifiers, api_base_url: str = STRING_API_BASE_DEFAULT):
    url = f"{api_base_url}/json/get_string_ids"
```

Thread through `determine_species()`:
```python
@staticmethod
def determine_species(df, nr=10, api_base_url: str = STRING_API_BASE_DEFAULT) -> int:
```

Thread through `get_species_taxon()`:
```python
def get_species_taxon(zip_path, df_list, nr=10, api_base_url: str = STRING_API_BASE_DEFAULT) -> int:
```

Import the default constant:
```python
from string_gsea.gsea_config import STRING_API_BASE_DEFAULT
```

### Step 0.5.6 — Update call site in `string_gsea_run.py`

**File:** `src/string_gsea/scripts/string_gsea_run.py`

Line 54, pass config's URL:
```python
species = get_species_taxon(zip_path, dataframes, api_base_url=config.api_base_url)
```

### Step 0.5.7 — Update `test_config_integration.py` (soon `test_config.py`)

Update the mocked API URL in test to use the constant if referenced.

### Step 0.5.8 — Verify

```bash
uv run pytest                              # unit tests pass
uv run pytest -m integration               # integration tests still work
uv run pytest -k test_gsea_session         # YAML round-trip includes api_base_url
```

---

## Files Modified (summary)

| File | Phase | Change |
|------|-------|--------|
| `pyproject.toml` | 0.1 | Add pytest markers config |
| `.gitignore` | 0.4 | Add `tests/.snakemake/` |
| `tests/test_get_species.py` | 0.2 | Add `@pytest.mark.integration` |
| `tests/test_gsea_results.py` | 0.2 | Replace skipif with `@pytest.mark.integration` |
| `tests/test_diff_xlsx_integration.py` | 0.3 | Rename to `test_diff_xlsx.py` |
| `tests/test_config_integration.py` | 0.3 | Rename to `test_config.py` |
| `src/string_gsea/gsea_config.py` | 0.5.1 | Add `api_base_url` field + constant |
| `src/string_gsea/gsea_session.py` | 0.5.2 | Replace class var with property |
| `src/string_gsea/string_gsea_builder.py` | 0.5.3 | Use config URL |
| `src/string_gsea/scripts/string_ora_run.py` | 0.5.4 | Thread config URL |
| `src/string_gsea/get_species.py` | 0.5.5 | Add `api_base_url` param |
| `src/string_gsea/scripts/string_gsea_run.py` | 0.5.6 | Pass config URL to `get_species_taxon` |
