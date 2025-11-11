"""
Getting Started with STRING-GSEA - Marimo Notebook

A simple introduction to using the STRING-GSEA Python package interactively.
"""

import marimo

__generated_with = "0.9.14"
app = marimo.App(width="medium")


@app.cell
def __():
    import marimo as mo
    return mo,


@app.cell
def __(mo):
    mo.md(
        """
        # Getting Started with STRING-GSEA
        
        Welcome to the **STRING-GSEA** interactive notebook! This notebook demonstrates 
        how to use the STRING-GSEA Python package for gene set enrichment analysis.
        
        ## What is STRING-GSEA?
        
        STRING-GSEA is a Python toolkit that:
        - Submits gene set enrichment analyses to STRING-DB API
        - Polls for completion automatically
        - Fetches and exports results (tables, graphs, links)
        - Generates interactive HTML reports
        """
    )
    return


@app.cell
def __(mo):
    mo.md(
        """
        ## Installation Check
        
        Let's verify that the required packages are installed:
        """
    )
    return


@app.cell
def __(mo):
    try:
        from string_gsea import __about__
        import polars as pl
        import plotly
        
        mo.md(f"""
        ✅ **string-gsea** version: `{__about__.__version__}`
        
        ✅ **polars** version: `{pl.__version__}`
        
        ✅ **plotly** version: `{plotly.__version__}`
        
        All required packages are installed!
        """)
    except ImportError as e:
        mo.md(f"""
        ❌ **Import Error:** {e}
        
        Please install the package with:
        ```bash
        uv pip install -e ".[notebooks]"
        ```
        """)
    return pl, plotly


@app.cell
def __(mo):
    mo.md(
        """
        ## Quick Example: Load GSEA Results
        
        Let's load some example GSEA results from the test data:
        """
    )
    return


@app.cell
def __():
    from pathlib import Path
    
    # Path to test data
    test_data_dir = Path("../../tests/data/WU_12345_GSEA")
    
    # Check if directory exists
    if test_data_dir.exists():
        excel_files = list(test_data_dir.glob("*.xlsx"))
        print(f"Found {len(excel_files)} Excel files")
    else:
        excel_files = []
        print(f"Test directory not found: {test_data_dir}")
    return Path, excel_files, test_data_dir


@app.cell
def __(excel_files, mo, pl):
    if excel_files:
        # Load the first Excel file
        df = pl.read_excel(excel_files[0])
        
        mo.md(f"""
        ### Loaded: `{excel_files[0].name}`
        
        **Shape:** {df.shape[0]} rows × {df.shape[1]} columns
        
        **First few rows:**
        """)
        
        mo.ui.table(df.head().to_pandas())
    else:
        df = None
        mo.md("⚠️ No test data found. Try running `string_gsea_run` first!")
    return df,


@app.cell
def __(mo):
    mo.md(
        """
        ## Interactive Controls
        
        Marimo notebooks are reactive! Changes to inputs automatically update downstream cells.
        
        Try adjusting the slider below:
        """
    )
    return


@app.cell
def __(mo):
    fdr_threshold = mo.ui.slider(
        start=0.01,
        stop=0.5,
        step=0.01,
        value=0.05,
        label="FDR Threshold",
        show_value=True
    )
    fdr_threshold
    return fdr_threshold,


@app.cell
def __(fdr_threshold, mo):
    mo.md(f"""
    You selected an FDR threshold of **{fdr_threshold.value}**
    
    This would filter pathways with FDR (False Discovery Rate) ≤ {fdr_threshold.value}
    """)
    return


@app.cell
def __(mo):
    mo.md(
        """
        ## Next Steps
        
        1. **Explore the GSEA Explorer notebook:** `gsea_exploration.py`
        2. **Check the documentation:** [STRING-GSEA README](../../README.md)
        3. **Run your own analysis:**
           ```bash
           string_gsea_run "path/to/data.zip" "workunit_id" "output_dir"
           ```
        4. **Generate reports:**
           ```bash
           string_gsea_render_reports "path/to/results"
           ```
        
        ## Resources
        
        - [Marimo Documentation](https://docs.marimo.io/)
        - [STRING-DB](https://string-db.org/)
        - [Polars Documentation](https://pola.rs/)
        """
    )
    return


if __name__ == "__main__":
    app.run()

