"""STRING-GSEA Index Report - Homepage with parameter summary and links."""

import marimo

__generated_with = "0.18.3"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _():
    import os
    import sys
    from pathlib import Path

    return Path, os, sys


@app.cell
def _(Path, os, sys):
    # Configuration: Try command-line args first, then environment variables
    if len(sys.argv) > 1:
        from string_gsea.marimo_reports._report_utils import ReportConfig

        config = ReportConfig.from_args()
    else:
        # Try environment variables
        data_file = os.environ.get("GSEA_DATA_FILE", "")
        links_file = os.environ.get("GSEA_LINKS_FILE", "")

        if data_file and links_file:
            from string_gsea.marimo_reports._report_utils import ReportConfig

            config = ReportConfig.from_env()
        else:
            # Fallback for interactive development
            from string_gsea.marimo_reports._report_utils import ReportConfig

            # Use a default path for development
            _base = Path(__file__).parent.parent.parent.parent
            _test_data = _base / "tests" / "data" / "fixtures" / "dummy_output"
            _gsea_dir = _test_data / "WU_abcd_GSEA" / "from_rnk"

            config = ReportConfig(
                data_file=_gsea_dir / "WUabcd_string_gsea_results_long.xlsx",
                links_file=_gsea_dir / "links.txt",
                fdr_threshold=0.05,
                genes_mapped_threshold=10,
            )

    return (config,)


@app.cell
def _(config):
    from string_gsea.marimo_reports._report_utils import get_workunit_id, load_links

    workunit_id = get_workunit_id(config)
    links = load_links(config)
    num_contrasts = len(links)

    return links, num_contrasts, workunit_id


@app.cell
def _(config, mo, workunit_id):
    mo.md(
        f"""
# STRING-GSEA Results: {workunit_id}

This report visualizes Gene Set Enrichment Analysis (GSEA) results
from [STRING-DB](https://string-db.org/).

## Parameters

| Parameter | Value |
|-----------|-------|
| FDR threshold | {config.fdr_threshold} |
| Genes mapped threshold | {config.genes_mapped_threshold} |
| Data file | `{config.data_file.name}` |
| Links file | `{config.links_file.name}` |
"""
    )
    return


@app.cell
def _(links, mo):
    # Create a table of contrasts with links
    _rows = []
    for contrast, url in links.items():
        _rows.append(f"| {contrast} | [{url}]({url}) |")

    _table = "\n".join(_rows)

    mo.md(
        f"""
## STRING-DB GSEA Results Links

| Contrast | URL |
|----------|-----|
{_table}
"""
    )
    return


@app.cell
def _(mo, num_contrasts):
    # Visualization reports with conditional multiple contrasts link
    _reports = ["- [Per-Contrast Visualizations (networks.html)](networks.html)"]

    if num_contrasts > 1:
        _reports.append("- [Cross-Contrast Comparisons (multiple.html)](multiple.html)")

    _reports_md = "\n".join(_reports)

    mo.md(
        f"""
## Visualization Reports

{_reports_md}
"""
    )
    return


if __name__ == "__main__":
    app.run()
