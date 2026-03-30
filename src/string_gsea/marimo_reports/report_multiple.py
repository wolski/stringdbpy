"""STRING-GSEA Multiple Contrasts Report - Cross-contrast comparisons (static)."""

import marimo

__generated_with = "0.18.3"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _():
    import io
    import os
    import sys
    import warnings
    from pathlib import Path

    import matplotlib.pyplot as plt
    import polars as pl

    # Suppress chained assignment warnings
    warnings.filterwarnings(
        "ignore",
        message="A value is trying to be set on a copy of a DataFrame",
        category=FutureWarning,
    )

    return Path, io, os, pl, plt, sys, warnings


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
def _():
    from string_gsea.dotplot_endrichment import dotplot_enrichment
    from string_gsea.gsea_plotting import make_upset_contrasts_terms
    from string_gsea.marimo_reports._report_utils import get_workunit_id, load_data
    from string_gsea.network import (
        add_gene_ratio,
        explode_protein_columns,
        filter_by_FDR,
        select_top_terms,
        summarize_terms,
    )

    return (
        add_gene_ratio,
        dotplot_enrichment,
        explode_protein_columns,
        filter_by_FDR,
        get_workunit_id,
        load_data,
        make_upset_contrasts_terms,
        select_top_terms,
        summarize_terms,
    )


@app.cell
def _(
    add_gene_ratio,
    config,
    explode_protein_columns,
    filter_by_FDR,
    get_workunit_id,
    load_data,
    select_top_terms,
    summarize_terms,
):
    # Load and process data
    df = load_data(config)
    workunit_id = get_workunit_id(config)

    df = filter_by_FDR(df, config.fdr_threshold, config.genes_mapped_threshold)
    df = select_top_terms(df, max_terms=config.max_terms)
    df = add_gene_ratio(df)
    xd = explode_protein_columns(df)
    xd = summarize_terms(xd)

    # Get unique categories
    category_list = sorted(xd.get_column("category").unique().to_list())

    return category_list, df, workunit_id, xd


@app.cell
def _(mo):
    # Backlink to index
    mo.md("[← Back to Index](index.html)")
    return


@app.cell
def _(config, mo, workunit_id):
    mo.md(
        f"""
# Cross-Contrast GSEA Comparisons: {workunit_id}

Comparing GSEA results across multiple contrasts with
FDR threshold = {config.fdr_threshold} and
genes mapped threshold = {config.genes_mapped_threshold}.

For each gene set category, we show:

1. **Upset plot**: Overlap of terms across contrasts
2. **Dotplot**: Enrichment across all contrasts
"""
    )
    return


@app.cell
def _(
    category_list,
    dotplot_enrichment,
    io,
    make_upset_contrasts_terms,
    mo,
    pl,
    plt,
    xd,
):
    def _render_category_plots(category):
        """Render upset plot and dotplot for a category."""
        plots = []

        # Filter data for this category
        xdf = xd.filter(pl.col("category") == category)

        if xdf.is_empty():
            return mo.md("_No data for this category_")

        # Check if we have at least 2 contrasts
        if xdf["contrast"].n_unique() < 2:
            return mo.md(
                "_Need at least 2 contrasts for cross-contrast comparison. "
                "Only 1 contrast found._"
            )

        # Upset plot
        try:
            make_upset_contrasts_terms(xd, category=category)
            _buf = io.BytesIO()
            plt.savefig(_buf, format="png", dpi=100, bbox_inches="tight")
            plt.close()
            _buf.seek(0)
            plots.append(mo.md("### Upset Plot"))
            plots.append(mo.image(_buf.read()))
            plots.append(
                mo.md(
                    "_Upset plot showing overlap of terms across contrasts. "
                    "Left bars show contrast size, top bars show unique/shared terms._"
                )
            )
        except Exception as e:
            plots.append(mo.md(f"_Upset plot error: {e}_"))

        # Dotplot
        try:
            # Prepare data for dotplot
            xd_smart = xd.drop(
                [col for col in xd.columns if col.startswith("protein")]
            ).unique()
            xd_smart = xd_smart.filter(pl.col("category") == category)

            dotplot_enrichment(xd_smart)
            _buf = io.BytesIO()
            plt.savefig(_buf, format="png", dpi=100, bbox_inches="tight")
            plt.close()
            _buf.seek(0)
            plots.append(mo.md("### Enrichment Dotplot"))
            plots.append(mo.image(_buf.read()))
            plots.append(
                mo.md(
                    "_Dotplot showing enrichment across contrasts. "
                    "Size=gene ratio, color=-log10(FDR), border=direction._"
                )
            )
        except Exception as e:
            plots.append(mo.md(f"_Dotplot error: {e}_"))

        return mo.vstack(plots)

    # Build accordion sections for each category
    category_sections = {}
    for _category in category_list:
        category_sections[_category] = _render_category_plots(_category)

    mo.accordion(category_sections)

    return (category_sections,)


if __name__ == "__main__":
    app.run()
