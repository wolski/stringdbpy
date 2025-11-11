"""
GSEA Results Explorer - Marimo Notebook

An interactive notebook for exploring STRING-GSEA enrichment results.
"""

import marimo

__generated_with = "0.16.5"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import polars as pl
    from pathlib import Path
    import plotly.express as px
    import plotly.graph_objects as go
    return (mo,)


@app.cell
def _(mo):
    mo.md(
        """
    # STRING-GSEA Results Explorer

    This interactive notebook allows you to explore Gene Set Enrichment Analysis (GSEA) 
    results from STRING-DB.

    ## 1. Load Your Data

    Upload or specify the path to your GSEA results Excel file:
    """
    )
    return


@app.cell
def _(mo):
    mo.md(
        """
    ---

    ## 4. Export Options

    You can export the filtered results for further analysis.
    """
    )
    return


if __name__ == "__main__":
    app.run()
