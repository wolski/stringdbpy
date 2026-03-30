"""STRING-GSEA Networks Report - Per-contrast visualizations (static)."""

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
    from pathlib import Path

    import matplotlib.pyplot as plt
    import polars as pl

    return Path, io, os, pl, plt, sys


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
    from string_gsea.cluster_genesets import (
        convert_to_binary,
        make_nested_dict,
        plot_term_distance_heatmap,
    )
    from string_gsea.gsea_plotting import make_upset, plot_term_ridges
    from string_gsea.marimo_reports._report_utils import (
        get_workunit_id,
        load_data,
        load_links,
    )
    from string_gsea.network import (
        add_gene_ratio,
        explode_protein_columns,
        filter_by_FDR,
        select_top_terms,
        summarize_terms,
    )
    from string_gsea.TermNetworkPlotter import TermNetworkPlotter

    return (
        TermNetworkPlotter,
        add_gene_ratio,
        convert_to_binary,
        explode_protein_columns,
        filter_by_FDR,
        get_workunit_id,
        load_data,
        load_links,
        make_nested_dict,
        make_upset,
        plot_term_distance_heatmap,
        plot_term_ridges,
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
    load_links,
    make_nested_dict,
    select_top_terms,
    summarize_terms,
):
    # Load and process data
    df = load_data(config)
    links = load_links(config)
    workunit_id = get_workunit_id(config)

    df = filter_by_FDR(df, config.fdr_threshold, config.genes_mapped_threshold)
    df = select_top_terms(df, max_terms=config.max_terms)
    df = add_gene_ratio(df)
    xd = explode_protein_columns(df)
    xd = summarize_terms(xd)
    nested = make_nested_dict(xd)

    # Get unique contrasts and categories
    contrast_list = sorted(xd.get_column("contrast").unique().to_list())
    category_list = sorted(xd.get_column("category").unique().to_list())

    return category_list, contrast_list, df, links, nested, workunit_id, xd


@app.cell
def _(mo):
    # Backlink to index
    mo.md("[← Back to Index](index.html)")
    return


@app.cell
def _(config, mo, workunit_id):
    mo.md(
        f"""
# Per-Contrast GSEA Visualizations: {workunit_id}

Visualizing STRING-DB GSEA results with FDR threshold = {config.fdr_threshold}
and genes mapped threshold = {config.genes_mapped_threshold}.

For each contrast and gene set category, we show:

1. **Gene Sets Table**: Summary of enriched terms
2. **Ridgeline plot**: Protein value distributions per term
3. **Upset plot**: Term/protein overlaps
4. **Heatmap**: Term and protein distances using Jaccard index
5. **Network plot**: Term-term relationships
"""
    )
    return


@app.cell
def _(
    TermNetworkPlotter,
    category_list,
    contrast_list,
    convert_to_binary,
    io,
    links,
    make_upset,
    mo,
    nested,
    pl,
    plot_term_distance_heatmap,
    plot_term_ridges,
    plt,
    xd,
):
    def _render_category_plots(contrast, category):
        """Render all plots for a contrast/category combination."""
        plots = []

        # Filter data
        xdf = xd.filter(
            (pl.col("contrast") == contrast) & (pl.col("category") == category)
        )

        if xdf.is_empty():
            return mo.md("_No data for this category_")

        # 1. Gene Sets Table (first)
        try:
            # Get unique gene sets with relevant columns
            _cols_to_show = ["term", "description", "fdr", "gene_ratio", "direction"]
            _available_cols = [c for c in _cols_to_show if c in xdf.columns]
            if _available_cols:
                gene_sets_df = xdf.select(_available_cols).unique().sort("fdr")
                plots.append(mo.md("#### Gene Sets"))
                plots.append(mo.ui.table(gene_sets_df.to_pandas(), selection=None))
        except Exception as e:
            plots.append(mo.md(f"_Gene sets table error: {e}_"))

        # 2. Ridgeline plot (second - shows distribution)
        try:
            plot_term_ridges(xdf)
            _buf = io.BytesIO()
            plt.savefig(_buf, format="png", dpi=100, bbox_inches="tight")
            plt.close()
            _buf.seek(0)
            plots.append(mo.md("#### Ridgeline Plot"))
            plots.append(mo.image(_buf.read()))
            plots.append(
                mo.md(
                    "_Distribution of protein input values per term. "
                    "Border=direction, bar height=-log10(FDR), color=gene ratio._"
                )
            )
        except Exception as e:
            plots.append(mo.md(f"_Ridgeline plot error: {e}_"))

        # 3. Upset plot (third - shows overlaps)
        _one_nested = nested.get(contrast, {}).get(category)
        if _one_nested is not None:
            _num_terms = len(_one_nested.columns) - 1
            if _num_terms >= 2 and len(_one_nested.columns) >= 3:
                try:
                    _binary = convert_to_binary(_one_nested, to_boolean=True)
                    _upset = make_upset(_binary, xdf)
                    _upset.plot()
                    _buf = io.BytesIO()
                    plt.savefig(_buf, format="png", dpi=100, bbox_inches="tight")
                    plt.close()
                    _buf.seek(0)
                    plots.append(mo.md("#### Upset Plot"))
                    plots.append(mo.image(_buf.read()))
                    plots.append(
                        mo.md(
                            "_Overlap of terms and proteins. "
                            "Left=term size, top=unique/shared proteins._"
                        )
                    )
                except Exception as e:
                    plots.append(mo.md(f"_Upset plot error: {e}_"))

            # 4. Heatmap/Clustering (fourth)
            if _num_terms >= 2:
                try:
                    _g = plot_term_distance_heatmap(_one_nested)
                    plt.setp(_g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
                    _buf = io.BytesIO()
                    plt.savefig(_buf, format="png", dpi=100, bbox_inches="tight")
                    plt.close()
                    _buf.seek(0)
                    plots.append(mo.md("#### Term-Protein Heatmap"))
                    plots.append(mo.image(_buf.read()))
                    plots.append(
                        mo.md(
                            "_Term and protein distances using Jaccard index. "
                            "Color=protein input value._"
                        )
                    )
                except Exception as e:
                    plots.append(mo.md(f"_Heatmap error: {e}_"))
            else:
                plots.append(
                    mo.md(
                        f"_Heatmap/Upset skipped: Only {_num_terms} term(s). Need at least 2._"
                    )
                )

        # 5. Network plot (last)
        try:
            _fig, _ax = plt.subplots(figsize=(10, 8))
            TermNetworkPlotter.plot_network(
                xd, category=category, contrast=contrast, thresh=1
            )
            _buf = io.BytesIO()
            plt.savefig(_buf, format="png", dpi=100, bbox_inches="tight")
            plt.close()
            _buf.seek(0)
            plots.append(mo.md("#### Network Plot"))
            plots.append(mo.image(_buf.read()))
            _legend = TermNetworkPlotter.get_figure_legend(category=category, thresh=1)
            plots.append(mo.md(_legend))
        except Exception as e:
            plots.append(mo.md(f"_Network plot error: {e}_"))

        return mo.vstack(plots)

    # Build accordion sections for each contrast
    contrast_sections = {}
    for _contrast in contrast_list:
        _url = links.get(_contrast.removesuffix("_results.tsv"), "")
        _link_md = f"[View in STRING-DB]({_url})" if _url else ""

        # Build tabs for each category within this contrast
        category_tabs = {}
        for _category in category_list:
            category_tabs[_category] = _render_category_plots(_contrast, _category)

        contrast_sections[_contrast] = mo.vstack(
            [mo.md(_link_md) if _link_md else mo.md(""), mo.ui.tabs(category_tabs)]
        )

    mo.accordion(contrast_sections)

    return (contrast_sections,)


if __name__ == "__main__":
    app.run()
