# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**STRING-GSEA** is a Python toolkit for submitting gene set enrichment analyses (GSEA) to the STRING-DB API, polling for completion, and fetching/exporting results. The project uses a builder pattern for clean separation between job submission and result handling.

Key dependencies: Python 3.13+, UV package manager, Quarto (for report generation)

## Development Commands

### Installation and Setup

```bash
# Install in editable mode (recommended for development)
uv pip install -e .

# Install with test dependencies
uv pip install -e ".[test]"

# Install with notebook dependencies (Marimo)
uv pip install -e ".[notebooks]"
```

### Configuration

```bash
# Create configuration file at $HOME/.config/string_gsea/config.toml
string_gsea_write_config --help
```

### Testing

```bash
# Run all tests with nox (uses uv backend)
nox -s test

# Run tests with pytest directly
pytest --durations=50 tests

# Run specific test file
pytest tests/test_gsea_session.py

# Run with additional pytest arguments via nox
nox -s test -- -v -k "test_name"
```

### Running Analysis

```bash
# Run GSEA with XLSX input (default)
string_gsea_run "path/to/data.zip" "workunit_id" "output_dir"

# Run GSEA with RNK files
string_gsea_run "path/to/data.zip" "workunit_id" "output_dir" --from-rnk

# Specify analysis type for XLSX input
string_gsea_run "path/to/data.zip" "workunit_id" "output_dir" --which "pep_1"
```

Valid analysis types: `pep_1`, `pep_1_no_imputed`, `pep_2`, `pep_2_no_imputed`

### Report Generation

```bash
# Render Quarto reports to output directory
string_gsea_render_reports "path/to/results" "path/to/output"

# Render reports into input directory
string_gsea_render_reports "path/to/results"

# Custom thresholds
string_gsea_render_reports "path/to/results" --FDR-threshold 0.01 --genes-mapped-threshold 5
```

### Interactive Notebooks

```bash
# Run Marimo notebooks interactively
marimo run docs/marimo_notebooks/getting_started.py
marimo run docs/marimo_notebooks/gsea_exploration.py

# Edit notebooks
marimo edit docs/marimo_notebooks/gsea_exploration.py
```

### Running Internal Scripts

```bash
# Run internal validation scripts using nox
nox -s run-internal-scripts

# Run specific modules
nox -s run-internal-scripts -- string_gsea.string_gsea_builder
```

## Architecture

### Core Components

The codebase follows a builder pattern with clear separation of concerns:

1. **StringGSEABuilder** (`string_gsea_builder.py`)
   - Orchestrates the entire GSEA workflow
   - Submits rank data to STRING-DB API
   - Polls for job completion with configurable timeouts
   - Builds StringGSEAResults objects
   - Entry point: `StringGSEABuilder(rank_dataframes, config, workunit_id, species, base_path)`

2. **GSEASession** (`gsea_session.py`)
   - Data class for session state management
   - Handles YAML serialization/deserialization of session data
   - Stores job IDs, results, and configuration
   - Keys for res_job_id and res_data are tuples: `(outer_key, inner_key)`
   - Methods: `to_yaml()`, `from_yaml()`

3. **StringGSEAResults** (`string_gsea_results.py`)
   - Handles downloading and writing of results from STRING-DB
   - Downloads TSV files, PNG graphs, and generates links
   - Creates output directory structure
   - Methods: `write_links()`, `write_gsea_tsv()`, `write_gsea_graphs()`, `zip_folder()`

4. **GSEAResultProcessor** (`gsea_result_processor.py`)
   - Post-processes TSV results into Excel reports
   - Creates three Excel formats: long, pivoted, and merged
   - Uses pyexcelerate for efficient Excel writing
   - Static method: `result_to_xlsx(tsv_dir, workunit_id)`

### Input Processing

5. **DiffXLSX** (`ranks_from_dea_xlsx.py`)
   - Extracts rank data from differential expression Excel files in ZIP archives
   - Filters by peptide counts and imputation status
   - Returns dict with tuple keys: `{(analysis_type, contrast): DataFrame}`
   - Analysis types: `pep_1`, `pep_1_no_imputed`, `pep_2`, `pep_2_no_imputed`

6. **GetTaxonID / OxFieldsZip** (`get_species.py`)
   - Determines species from protein identifiers
   - Uses STRING API to fetch NCBI Taxon IDs
   - Contains embedded species mapping data (ZIP files in `src/string_gsea/data/mappings/`)
   - Function: `get_species_taxon(zip_path, dataframes)` returns taxon ID (e.g., 9606 for human)

### Configuration

7. **GSEAConfig** (`gsea_config.py`)
   - Data class for GSEA configuration
   - Fields: `api_key`, `fdr`, `caller_identity`, `ge_enrichment_rank_direction`
   - Reads from `$HOME/.config/string_gsea/config.toml`
   - Function: `get_configuration()` loads config

### Visualization and Reporting

8. **Quarto Reports** (`src/string_gsea/docs/`)
   - Template-based reporting using Quarto
   - Profiles: `multiple` (default), `single`
   - Templates: `VisualizeMultipleContrastsGSEA.qmd`, `VisualizeNetworks.qmd`
   - Configuration: `_quarto.yml`
   - Rendered via `render_reports.py` which calls `quarto render` with parameters

9. **Network Visualization** (`network.py`, `TermNetworkBuilder.py`, `TermNetworkPlotter.py`)
   - Creates network graphs from GSEA results
   - Uses networkx for graph structure
   - Supports Cytoscape and pyvis visualization

10. **Plotting** (`gsea_plotting.py`, `dotplot_endrichment.py`)
    - Matplotlib/seaborn visualizations
    - Dotplots for enrichment results

### Scripts (Entry Points)

Located in `src/string_gsea/scripts/`:
- `string_gsea_run.py`: Main analysis workflow (uses cyclopts for CLI)
- `render_reports.py`: Quarto report rendering
- `write_config.py`: Configuration file creation

All scripts are registered in `pyproject.toml` under `[project.scripts]`

## Data Flow

1. **Input**: ZIP file containing either RNK files or XLSX with differential expression data
2. **Species Detection**: Automatic taxon ID detection from identifiers
3. **Rank Extraction**: Parse input into rank DataFrames (dict with tuple keys)
4. **Job Submission**: Submit to STRING-DB API, receive job IDs
5. **Polling**: Wait for completion (configurable timeout, default 3600s)
6. **Results Download**: Fetch TSV, PNG, and links from STRING-DB
7. **Post-Processing**: Generate Excel reports (long, pivoted, merged)
8. **Session Serialization**: Save YAML session file for reproducibility
9. **Report Rendering**: Optional Quarto HTML reports with visualizations

## Output Structure

After running `string_gsea_run`:
```
output_directory/
├── WU_{workunit_id}_GSEA/
│   ├── {contrast_name}/
│   │   ├── {inner_key}_results.tsv
│   │   ├── {inner_key}_results.png
│   │   ├── links.txt
│   │   └── *.rnk
│   ├── gsea_session.yml
│   ├── WU_{workunit_id}_string_gsea_results_long.xlsx
│   ├── WU_{workunit_id}_string_gsea_results_pivoted.xlsx
│   └── WU_{workunit_id}_string_gsea_results_merged.xlsx
└── WU_{workunit_id}_GSEA.zip
```

After running `render_reports`:
```
output_directory/rendered_reports/
├── index.html
├── VisualizeNetworks.html
└── EnrichmentResults.html
```

## Important Implementation Details

### Nested Dictionary Keys
Throughout the codebase, results are keyed by tuples `(outer_key, inner_key)`:
- **outer_key**: Analysis type (e.g., "pep_1", "pep_2_no_imputed") or contrast group
- **inner_key**: Specific contrast name

When serializing to YAML, these are converted to strings with "~" separator: `"outer~inner"`

### API Configuration
STRING-DB API base: `https://version-12-0.string-db.org/api`
Endpoints:
- Submit: `/json/valuesranks_enrichment_submit`
- Status: `/json/valuesranks_enrichment_status`
- Get IDs: `/json/get_string_ids`

Required parameters: `species`, `caller_identity`, `api_key`, `ge_fdr`, `ge_enrichment_rank_direction`

### Polling Behavior
Default polling: 10-second intervals, 3600-second max timeout
Failure statuses: `'nothing found'`, `'unknown organism'`
Success status: `'success'`

### Data Processing with Polars
The project uses Polars DataFrames exclusively (not pandas). When adding features:
- Use `pl.read_csv()`, `pl.read_excel()` for reading
- Use `.write_csv()`, `.select()`, `.filter()` for transformations
- Results include computed column: `directionNR` (1 for "top", -1 for "bottom", 0 otherwise)

## Testing Notes

Test data located in `tests/data/`:
- `DE_mouse_fasta_xlsx.zip`: Mouse differential expression example
- `2848501.zip`: Human RNK file example

Tests use pytest with fixtures. When writing tests:
- Mock external API calls to STRING-DB
- Use temporary directories for file outputs
- Test YAML serialization round-trips
