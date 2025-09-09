# stringdbpy

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.13+](https://img.shields.io/badge/python-3.13+-blue.svg)](https://www.python.org/downloads/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

**STRING‑GSEA** is a Python toolkit to submit gene set enrichment analyses (GSEA) to the STRING‑DB API, poll for completion, and fetch/export results (tables, graphs, links) with minimal boilerplate.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Command‑Line Interface](#command-line-interface)
- [Python API Usage](#python-api-usage)
- [Configuration](#configuration)
- [Session Workflow](#session-workflow)
- [Examples](#examples)
- [Development](#development)
- [License](#license)

---

## Features

- **Builder pattern** for clean separation between job submission and result handling via `StringGSEABuilder`
- **YAML‑backed session** to serialize/deserialize STRING‑DB session using `GSEASession`
- **Automatic polling** of STRING‑DB until results are ready with configurable timeouts
- **Multiple input formats** support for both RNK files and XLSX files with differential expression data
- **Export utilities** for rank files, result TSVs, enrichment graphs, and links via `StringGSEAResults`
- **Excel report generation** with pivoted and merged formats via `GSEAResultProcessor`
- **Zip archiving** of full result folders
- **Command-line scripts** for one‑step runs (`string_gsea_run`) and report rendering (`render_reports`)
- **Quarto report generation** with interactive visualizations and parameterized templates
- **Species detection** from input files using `get_species_taxon`
- **Configurable parameters** for FDR thresholds, enrichment direction, and API settings

## Installation

### Prerequisites

- Python 3.13 or higher
- [UV](https://docs.astral.sh/uv/) package manager (recommended) or pip
- [Quarto](https://quarto.org/) (for report generation)

### Installation Options

#### Option 1: Development Installation (Recommended)

```bash
# Clone the repository
git clone https://github.com/wolski/stringdbpy.git
cd stringdbpy

# Install in editable mode with UV
uv pip install -e .

# Or with pip
pip install -e .
```

#### Option 2: Install from PyPI (when available)

```bash
uv pip install string-gsea
# or
pip install string-gsea
```

### Post-Installation Setup

1. **Create configuration file:**
   ```bash
   string_gsea_write_config --help
   ```
   This creates a configuration file at `$HOME/.config/string_gsea/config.toml` containing your API key, FDR threshold, and caller identity.

2. **Verify installation:**
   ```bash
   string_gsea_run --help
   render_reports --help
   ```

### Dependencies

The package automatically installs all required dependencies including:
- `cyclopts` for command-line interface
- `polars` for data processing
- `requests` for API communication
- `loguru` for logging
- `pyyaml` for configuration
- `quarto` for report generation
- And other scientific computing libraries

## Command‑Line Interface

The package provides three main command-line tools for different aspects of the STRING-GSEA workflow.

### 1. Configuration Setup

First, set up your configuration file:

```bash
string_gsea_write_config --help
```

This creates a configuration file at `$HOME/.config/string_gsea/config.toml` containing:
- `api_key`: Your STRING-DB API key
- `fdr`: False discovery rate threshold (default: 0.25)
- `caller_identity`: Unique identifier for API calls
- `ge_enrichment_rank_direction`: Direction for enrichment ranking (1 or -1)

### 2. Running GSEA Analysis

The main analysis command:

```bash
string_gsea_run --help
```

#### Basic Usage

```bash
# Run analysis with XLSX input (default)
string_gsea_run "path/to/data.zip" "workunit_id" "output_directory"

# Run analysis with RNK files
string_gsea_run "path/to/data.zip" "workunit_id" "output_directory" --from-rnk

# Specify analysis type for XLSX input
string_gsea_run "path/to/data.zip" "workunit_id" "output_directory" --no-from-rnk --which "pep_1"
```

#### Parameters

- `zip_path`: Path to input ZIP file containing rank data
- `workunit_id`: Unique identifier for this analysis run
- `out_dir`: Output directory (default: current directory)
- `--from-rnk`: Use RNK files instead of XLSX (default: False)
- `--which`: Analysis type for XLSX input (default: "pep_2_no_imputed")

#### Examples

```bash
# Mouse data analysis
string_gsea_run "./tests/data/DE_mouse_fasta_xlsx.zip" "12345" "./results/mouse_analysis"

# Human data with RNK files
string_gsea_run "./tests/data/2848501.zip" "abcd" "./results/human_analysis" --from-rnk

# Custom analysis type
string_gsea_run "./data/experiment.zip" "exp001" "./results" --which "pep_1"
```

### 3. Generating Reports

Create interactive HTML reports from analysis results:

```bash
string_gsea_run --help
string_gsea_render_reports --help
```

#### Basic Usage

```bash
# Render reports to a specific output directory
render_reports "path/to/results" "path/to/output"

# Render reports into the input directory
render_reports "path/to/results"

# Custom FDR threshold and gene mapping threshold
render_reports "path/to/results" --FDR-threshold 0.01 --genes-mapped-threshold 5
```

#### Parameters

- `data_dir`: Path to results directory (containing XLSX and links files)
- `output_dir`: Output directory for rendered reports (optional)
- `--FDR-threshold`: FDR threshold for filtering (default: 0.05)
- `--genes-mapped-threshold`: Minimum genes mapped threshold (default: 10)

#### Examples

```bash
# Standard report generation
render_reports "./results/mouse_analysis" "./reports/mouse_reports"

# Render into results directory
render_reports "./results/human_analysis"

# Custom thresholds
render_reports "./results/experiment" --FDR-threshold 0.01 --genes-mapped-threshold 15
```

### Output Structure

After running `string_gsea_run`, you'll find:

```
output_directory/
├── WU_workunit_id_GSEA/
│   ├── contrast1/
│   │   ├── results.tsv
│   │   ├── results.png
│   │   ├── links.txt
│   │   └── *.rnk
│   ├── contrast2/
│   │   └── ...
│   ├── gsea_session.yml
│   ├── WU_workunit_id_string_gsea_results_long.xlsx
│   ├── WU_workunit_id_string_gsea_results_pivoted.xlsx
│   └── WU_workunit_id_string_gsea_results_merged.xlsx
└── WU_workunit_id_GSEA.zip
```

After running `render_reports`, you'll find:

```
output_directory/
└── rendered_reports/
    ├── index.html
    ├── VisualizeNetworks.html
    ├── EnrichmentResults.html
    └── ...
```

## TODO

### High Priority

- **integration with bfabric**
- **R Integration**: Create R functions to convert STRING-DB GSEA results to ClusterProfiler R objects for seamless integration with the R bioinformatics ecosystem
- **scverse integration** Improve interoperability with scverse (e.g. decoupler)


### Medium Priority

- **Interactive Visualizations**: Port static visualizations from matplotlib/seaborn to Altair for enhanced interactivity in HTML reports

### Low Priority

- **Documentation**: Expand documentation with more examples and tutorials


## License

Distributed under the [MIT License](https://spdx.org/licenses/MIT.html).

