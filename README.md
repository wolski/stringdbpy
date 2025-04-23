# A373_STRING_GSEA

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

- **Builder pattern** for clean separation between job submission and result handling
- **YAML‑backed session** to serialize/deserialize in‑flight or completed runs
- **Automatic polling** of STRING‑DB until results are ready
- **Export utilities** for rank files, result TSVs, enrichment graphs, and links
- **Zip archiving** of full result folders
- **Console script** for one‑step runs from the command line

## Installation

```shell
# Clone the repo and install in editable mode
git clone https://github.com/yourorg/A373_STRING_GSEA.git
cd A373_STRING_GSEA
pip install -e .
```

## Quick Start

```python
from pathlib import Path
from string_gsea.string_gsea_builder import StringGSEABuilder
from string_gsea.gsea_utilities import get_rank_files

# Prepare inputs
zip_path = Path("/path/to/your/rank_files.zip")
workunit = "WU123"
config = {
  "api_key": "YOUR_API_KEY",
  "fdr": 0.25,
  "caller_identity": "your.tool.name",
  "ge_enrichment_rank_direction": -1
}

# Run builder & poll
builder = StringGSEABuilder(
  rank_dataframes=get_rank_files(zip_path),
  config_dict=config,
  workunit_id=workunit,
  species=9606,
  base_path=Path("./output")
)
results = builder.get_result()
builder.save_session()               # YAML

# Export outputs
results.write_links()
results.write_gsea_tsv()
results.write_gsea_graphs()
```

## Command‑Line Interface

Once installed, you can run from the shell:

```bash
string-gsea-builder \
  tests/data/2848501.zip \
  --workunit-id WU123 \
  --fdr 0.25 \
  --out-dir ./results
```

## Python API Usage

1. **Instantiate the builder**: supply rank DataFrames and config.  
2. **Write rank files** (optional): `builder.write_rank_files()`  
3. **Submit & poll**: `builder.submit().poll()` or `builder.get_result()`  
4. **Serialize session**: `builder.save_session()`  
5. **Build results object**: `builder.build_results()`  
6. **Export**: use `StringGSEAResults` methods to write links, TSVs, graphs, and zip archives.

## Configuration

Your `config_dict` must include:

- `api_key`: your STRING‑DB API key  
- `fdr`: false discovery rate cutoff  
- `caller_identity`: a unique caller name  
- `ge_enrichment_rank_direction`: `1` or `-1`

## Session Workflow

Outlines the lifecycle of a run:

1. **YAML session**: saved by `builder.save_session()`  
2. **Reload**: `StringGSEABuilder.load_session(path)` or `GSEASession.from_yaml()`  
3. **Continue**: if polling not done, `.poll()`, then export results.

## Examples

Additional examples and use cases will be added here.

## Development

1. Fork & clone this repo  
2. Create a virtual environment and install dependencies  
3. Run tests: `pytest`  
4. Lint: `flake8` / `black .`  

## License

Distributed under the [MIT License](https://spdx.org/licenses/MIT.html).

