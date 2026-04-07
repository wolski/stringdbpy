# stringdbpy

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.13+](https://img.shields.io/badge/python-3.13+-blue.svg)](https://www.python.org/downloads/)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)

**STRING-GSEA** is a Python toolkit to submit gene set enrichment analyses (GSEA) and over-representation analyses (ORA) to the STRING-DB API, poll for completion, and fetch/export results (tables, graphs, links) with minimal boilerplate.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Command-Line Interface](#command-line-interface)
- [Batch Workflow](#batch-workflow)
- [Docker](#docker)
- [Interactive Notebooks](#interactive-notebooks)
- [TODO](#todo)
- [License](#license)

---

## Features

- **Builder pattern** for clean separation between job submission and result handling via `StringGSEABuilder`
- **YAML-backed session** to serialize/deserialize STRING-DB session using `GSEASession`
- **Automatic polling** of STRING-DB until results are ready with configurable timeouts
- **Multiple input formats** support for both RNK files and XLSX files with differential expression data
- **ORA analysis** with custom background gene sets via `string_ora_run`
- **Export utilities** for rank files, result TSVs, enrichment graphs, and links via `StringGSEAResults`
- **Excel report generation** with pivoted and merged formats via `GSEAResultProcessor`
- **Quarto report generation** with interactive visualizations and parameterized templates
- **Snakemake batch workflow** for processing multiple datasets with a single command
- **Species detection** from input files using `get_species_taxon`
- **Configurable parameters** for FDR thresholds, enrichment direction, and API settings

## Installation

### Prerequisites

- Python 3.13 or higher
- [UV](https://docs.astral.sh/uv/) package manager (recommended) or pip
- [Quarto](https://quarto.org/) (for report generation)
- R with `stringGSEAplot` package (for Quarto reports)

### Install from GitHub

```bash
git clone https://github.com/wolski/stringdbpy.git
cd stringdbpy

# Install in editable mode with UV
uv pip install -e .

# Or with pip
pip install -e .
```

### Post-Installation Setup

Create a configuration file:

```bash
string_gsea_write_config
```

This creates `$HOME/.config/string_gsea/config.toml` containing:
- `api_key`: Your STRING-DB API key
- `fdr`: False discovery rate threshold (default: 0.25)
- `caller_identity`: Unique identifier for API calls
- `ge_enrichment_rank_direction`: Direction for enrichment ranking (1 or -1)

Verify installation:

```bash
string_gsea_run --help
string_ora_run --help
string_gsea_workflow --help
```

## Command-Line Interface

The package provides four command-line tools:

### `string_gsea_write_config` -- Configuration Setup

```bash
string_gsea_write_config
```

### `string_gsea_run` -- GSEA Analysis

```bash
# Run analysis with XLSX input (default)
string_gsea_run "path/to/data.zip" "workunit_id" "output_directory"

# Run analysis with RNK files
string_gsea_run "path/to/data.zip" "workunit_id" "output_directory" --from-rnk

# Specify analysis type for XLSX input
string_gsea_run "path/to/data.zip" "workunit_id" "output_directory" --which "pep_1"
```

Parameters:
- `zip_path`: Path to input ZIP file containing rank data
- `workunit_id`: Unique identifier for this analysis run
- `out_dir`: Output directory (default: current directory)
- `--from-rnk`: Use RNK files instead of XLSX (default: False)
- `--which`: Analysis type for XLSX input (default: "pep_2_no_imputed"). Valid: `pep_1`, `pep_1_no_imputed`, `pep_2`, `pep_2_no_imputed`

### `string_ora_run` -- Over-Representation Analysis

```bash
# Run ORA with significant genes, background, and FASTA for species detection
string_ora_run "significant.txt" "background.txt" "proteome.fasta" \
  --out-dir "./results" --workunit-id "ORA001"
```

Input files:
- `significant.txt`: One protein/gene ID per line (significant genes)
- `background.txt`: One protein/gene ID per line (all tested genes)
- `proteome.fasta`: FASTA file with OX= fields for species detection

### `string_gsea_workflow` -- Batch Snakemake Workflow

See [Batch Workflow](#batch-workflow) below.

### Output Structure

After running `string_gsea_run`:

```
output_directory/
├── WU_{workunit_id}_GSEA/
│   ├── {contrast_name}/
│   │   ├── results.tsv
│   │   ├── results.png
│   │   ├── links.txt
│   │   └── *.rnk
│   ├── gsea_session.yml
│   ├── WU_{workunit_id}_string_gsea_results_long.xlsx
│   ├── WU_{workunit_id}_string_gsea_results_pivoted.xlsx
│   └── WU_{workunit_id}_string_gsea_results_merged.xlsx
└── WU_{workunit_id}_GSEA.zip
```

## Batch Workflow

The `string_gsea_workflow` command runs GSEA, renders Quarto reports, and packages results for multiple datasets using Snakemake.

```bash
# Run all datasets
string_gsea_workflow --cores 4

# CI-tagged datasets only
string_gsea_workflow --target ci --cores 2

# Custom input/output directories
string_gsea_workflow --datasets-dir /path/to/datasets --output-dir /path/to/outputs

# Dry-run
string_gsea_workflow --dry-run

# Generate config
string_gsea_workflow config
```

The workflow expects datasets in `data/datasets/{name}/` with `input.zip` and `params.yml`. Results are written to `data/outputs/{name}/`.

## Docker

### Build

```bash
# Single-platform (current architecture)
docker buildx build -f docker/Dockerfile -t string-gsea:local --load .

# Multi-arch
docker buildx build -f docker/Dockerfile -t string-gsea:local --platform linux/amd64,linux/arm64 .
```

### Usage via wrapper script

The recommended way to run the container is via `docker/string_gsea_docker.sh`. It handles image pulling, `--init` for signal handling, `--user` so output files are owned by you, and mounts `$(pwd)` as `/work`. Prefers Podman if available, falls back to Docker.

```bash
# First time: generate config.toml (mounts ~/.config/string_gsea writable)
./docker/string_gsea_docker.sh config

# Run the pipeline from your project directory
# Expects data/datasets/{name}/input.zip + params.yml under cwd by default
cd /path/to/my/project
./docker/string_gsea_docker.sh --cores 4

# Or specify where your datasets and outputs live (relative to cwd)
./docker/string_gsea_docker.sh --datasets-dir my_datasets --output-dir my_results --cores 4

# Dry-run: show what would execute without running anything
./docker/string_gsea_docker.sh --dry-run

# Pin a specific image version
./docker/string_gsea_docker.sh --image-version 0.1.0 -- --cores 4
```

### Usage with plain docker

```bash
# Generate config (first time only)
docker run --rm -it --init \
  --user "$(id -u):$(id -g)" \
  -v ~/.config/string_gsea:/root/.config/string_gsea \
  -v "$(pwd):/work" -w /work \
  string-gsea:local config

# Run full pipeline (current directory mounted as /work)
docker run --rm -it --init \
  --user "$(id -u):$(id -g)" \
  -v "$(pwd):/work" -w /work \
  -v ~/.config/string_gsea:/root/.config/string_gsea:ro \
  string-gsea:local --cores 4

# Custom input/output directories
docker run --rm -it --init \
  --user "$(id -u):$(id -g)" \
  -v "$(pwd):/work" -w /work \
  -v ~/.config/string_gsea:/root/.config/string_gsea:ro \
  string-gsea:local --datasets-dir my_datasets --output-dir my_results --cores 4
```

The container expects the working directory to contain `data/datasets/{name}/input.zip` + `params.yml` (or use `--datasets-dir`/`--output-dir` to override). Results are written to `data/outputs/` by default.

## Interactive Notebooks

STRING-GSEA includes interactive [Marimo](https://marimo.io/) notebooks for exploratory analysis.

```bash
# Install notebook dependencies
uv pip install -e ".[notebooks]"

# Launch notebooks
marimo run docs/marimo_notebooks/getting_started.py
marimo run docs/marimo_notebooks/gsea_exploration.py

# Edit notebooks
marimo edit docs/marimo_notebooks/gsea_exploration.py
```

Available notebooks:
- **`getting_started.py`** -- Introduction to STRING-GSEA and Marimo basics
- **`gsea_exploration.py`** -- Interactive explorer for GSEA results with filtering and visualization

For detailed instructions, see the [Marimo Quick Start](docs/marimo_notebooks/QUICKSTART.md).

## TODO

### High Priority

- **Integration with bfabric**
- **R Integration**: Create R functions to convert STRING-DB GSEA results to ClusterProfiler R objects for seamless integration with the R bioinformatics ecosystem
- **scverse integration**: Improve interoperability with scverse (e.g. decoupler)

### Medium Priority

- **Interactive Visualizations**: Port static visualizations from matplotlib/seaborn to Altair for enhanced interactivity in HTML reports

### Low Priority

- **Documentation**: Expand documentation with more examples and tutorials

## License

Distributed under the [MIT License](https://spdx.org/licenses/MIT.html).
