# Marimo Notebooks Quick Start

## Install

```bash
uv pip install -e ".[notebooks]"
```

## Run

```bash
# Getting started tutorial
marimo edit docs/marimo_notebooks/getting_started.py

# GSEA results explorer
marimo edit docs/marimo_notebooks/gsea_exploration.py
```

## Key Commands

```bash
marimo edit notebook.py   # Edit mode (interactive development)
marimo run notebook.py    # Run mode (read-only web app)
marimo tutorial           # Interactive Marimo tutorial
```

## Example: Create New Notebook

```bash
marimo edit docs/marimo_notebooks/my_analysis.py
```

Then add cells like:

```python
import marimo as mo
import polars as pl

# Load data
df = pl.read_excel("results.xlsx")

# Interactive filter
fdr = mo.ui.slider(0.01, 0.5, value=0.05, label="FDR")
fdr

# Reactive filtering (auto-updates when slider changes!)
filtered = df.filter(pl.col('fdr') <= fdr.value)
mo.ui.table(filtered.to_pandas())
```

## Resources

- [Marimo Docs](https://docs.marimo.io/)
- [STRING-GSEA README](../../README.md)
