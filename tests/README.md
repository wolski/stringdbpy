# Tests

## Running Tests

From the project root:

```bash
make help              # show all targets
make test              # fast/unit tests only (default, no network)
make test-integration  # integration tests (requires STRING-DB access)
make test-all          # all tests (unit + integration)
make test-ci           # Snakemake CI datasets (small, single-contrast, requires STRING-DB)
make lint              # ruff lint check
make format            # ruff auto-format
make check             # lint + unit tests
```

Or directly with pytest:

```bash
uv run pytest tests                    # unit tests (integration skipped by default)
uv run pytest -m integration tests     # integration tests only
uv run pytest -m "" tests              # all tests
uv run pytest -k test_gsea_session     # specific test
```

## Test Types

### Unit tests (`make test`)

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

### Integration tests (`make test-integration`)

Require live STRING-DB API access. Marked with `@pytest.mark.integration`:

- `test_get_species.py::test_GetTaxonID_determine_species` — queries STRING API to determine species from protein identifiers
- `test_gsea_results.py::test_gsea_results_integration` — downloads results (TSV, PNG, links) from a saved session

### Snakemake CI tests (`make test-ci`)

End-to-end integration via Snakemake. Runs `string_gsea_run` + `string_gsea_render_marimo` on small datasets tagged `ci`. Must be run from `tests/` directory (the Makefile handles this).

```bash
cd tests && uv run snakemake -s Snakefile ci -j1      # CI datasets
cd tests && uv run snakemake -s Snakefile -j1          # all datasets
cd tests && uv run snakemake -s Snakefile run_mouse_xlsx  # single dataset
cd tests && uv run snakemake -s Snakefile -n           # dry-run
cd tests && uv run snakemake -s Snakefile help          # list datasets
```

## Directory Layout

```
tests/
├── README.md
├── conftest.py              # Shared pytest fixtures
├── Snakefile                # Snakemake workflow for batch integration tests
├── test_*.py                # Test files
└── data/
    ├── datasets/            # Input datasets (Snakemake-discoverable)
    │   ├── human_rnk_2848501/    # [ci] Human, RNK, 8.4 MB
    │   ├── human_rnk_2881688/    #       Human, RNK, 76 MB
    │   ├── human_rnk_2958180_single/ #   Human, RNK, 211 MB, single contrast
    │   ├── human_xlsx_20250813/  #       Human, XLSX, 44 MB
    │   ├── human_xlsx_20250908/  #       Human, XLSX, 39 MB
    │   ├── human_xlsx_primary_58167/ #   Human, XLSX, 80 MB
    │   ├── human_xlsx_secondary_58167/ # Human, XLSX, 82 MB
    │   ├── mouse_xlsx/           # [ci] Mouse, XLSX, 9.5 MB
    │   └── yeast_rnk/            # [ci] Yeast, RNK, 2.0 MB
    ├── fixtures/            # Test-specific fixtures (not for Snakemake)
    │   ├── dummy_output/    # Pre-generated XLSX for network tests
    │   ├── dummy_session/   # Saved session YAML + ZIP for results tests
    │   └── no_matching_zip/ # FASTA ZIP for species detection tests
    └── outputs/             # Snakemake output directory (gitignored)
```

Each dataset directory contains:
- `input.zip` — the input data (RNK files or XLSX with FASTA)
- `params.yml` — configuration: species, input type, workunit ID, tags, filtering

Datasets tagged `ci` are small (<10 MB) and used for fast integration testing.

## Outputs

- **pytest**: Results printed to terminal. No persistent output files (uses `tmp_path`).
- **Snakemake**: Output written to `tests/data/outputs/{dataset_name}/`. This directory is not checked into git. Clean with `cd tests && snakemake clean`.
