# Tests

## Running Tests

From the project root:

```bash
uv run pytest tests                    # unit tests (integration skipped by default)
uv run pytest -m integration tests     # integration tests only
uv run pytest -m "" tests              # all tests
uv run pytest -k test_gsea_session     # specific test
```

## Test Types

### Unit tests (default)

No network access. Run fast (~3s). These test:

- Config TOML read/write and validation (`test_config.py`)
- XLSX rank extraction (`test_diff_xlsx.py`)
- Species detection from FASTA OX= fields (`test_get_species.py`, local tests only)
- GSEA session YAML serialization round-trips (`test_gsea_session.py`)
- Utility functions (`test_gsea_utilities.py`)
- Network graph operations (`test_network.py`)
- ZIP/directory input handling (`test_prep_data_input.py`)
- Quarto report rendering (mocked subprocess) (`test_render_reports.py`)
- Taxon ID mapping from bundled data files (`test_taxon_utils.py`)

### Integration tests (`pytest -m integration`)

Require live STRING-DB API access. Marked with `@pytest.mark.integration`:

- `test_get_species.py::test_GetTaxonID_determine_species` — queries STRING API to determine species from protein identifiers
- `test_gsea_results.py::test_gsea_results_integration` — downloads results (TSV, PNG, links) from a saved session
- `test_workflow_integration.py::test_string_gsea_run` — runs `string_gsea_run` end-to-end on the `mouse_xlsx` dataset

## Directory Layout

```
tests/
├── README.md
├── conftest.py              # Shared pytest fixtures
├── test_*.py                # Test files
└── data/
    ├── datasets/            # Input datasets
    │   ├── human_rnk_2848501/
    │   ├── human_rnk_2881688/
    │   ├── human_rnk_2958180_single/
    │   ├── human_xlsx_20250813/
    │   ├── human_xlsx_20250908/
    │   ├── human_xlsx_primary_58167/
    │   ├── human_xlsx_secondary_58167/
    │   ├── mouse_xlsx/           # XLSX, 9.5 MB — used in integration tests
    │   └── yeast_rnk/
    ├── fixtures/            # Test-specific fixtures
    │   ├── dummy_output/    # Pre-generated XLSX for network tests
    │   ├── dummy_session/   # Saved session YAML + ZIP for results tests
    │   └── no_matching_zip/ # FASTA ZIP for species detection tests
    └── outputs/             # (gitignored)
```

Each dataset directory contains:
- `input.zip` — the input data (RNK files or XLSX with FASTA)
- `params.yml` — metadata (species, input type, workunit ID)

## Outputs

- **pytest**: Results printed to terminal. No persistent output files (uses `tmp_path`).
