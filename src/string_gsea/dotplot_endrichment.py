import polars as pl
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import altair as alt
import numpy as np

def prepare_data_for_plotting(pdf:pd.DataFrame, label_length:int=40):
    pdf['-log10FDR'] = -pdf['falseDiscoveryRate'].apply(lambda x: np.log10(x + 1e-10))

    direction_colors = {
        'bottom': 'blue',
        'top': 'yellow',
        'both ends': 'green'
    }
    pdf['borderColor'] = pdf['direction'].map(direction_colors).fillna('grey')

    # Size based on geneRatio
    pdf['circleSize'] = pdf['geneRatio'] * 500  # No normalization
    # Truncate termDescription for display
    # Create termLabel: "termID: truncated_description"
    def make_label(row):
        desc = row['termDescription']
        truncated = desc if len(desc) <= label_length else desc[:label_length - 1] + "..."
        return f"{row['termID']}: {truncated}"

    pdf['termLabel'] = pdf.apply(make_label, axis=1)
    return pdf, direction_colors

def plot_enrichment_scatter(ax, pdf) -> plt.scatter:
    # Use vectorized scatter for all points
    sc = ax.scatter(
        x=pdf['contrast'],
        y=pdf['termLabel'],
        s=pdf['circleSize'],
        c=pdf['-log10FDR'],
        cmap='coolwarm',
        vmin=pdf['-log10FDR'].min(),
        vmax=pdf['-log10FDR'].max(),
        edgecolors=pdf['borderColor'],
        linewidth=2
    )
    
    pdf['contrast_code'] =  pdf['contrast'].astype('category').cat.codes.astype(float)
    codes = pdf['contrast_code'].unique()
    ax.set_xlim(codes.min() - 0.5, codes.max() + 0.5)
    ax.set_xticks(codes)
    ax.set_xticklabels(pdf['contrast'].astype('category').cat.categories, rotation=45, ha='right')

    ax.margins(y=0)
    # 2) explicitly clamp the y‐limits to your term count
    n_terms = pdf['termLabel'].nunique()
    # if you've inverted the y‐axis to show first term on top:
    ax.set_ylim(n_terms - 0.5, -0.5)

    ax.set_xlabel('Contrast')
    ax.set_ylabel('Gene Set Terms')
    return sc

def add_custom_legends(fig, ax, direction_colors):
    # Border legend (Direction)
    border_handles = [
        mlines.Line2D([], [], marker='o', linestyle='None', markersize=10,
                      markerfacecolor='white', markeredgecolor=color, label=label)
        for label, color in direction_colors.items()
    ]
    legend1 = fig.legend(
        handles=border_handles,
        title='Direction',
        loc='center right',
        bbox_to_anchor=(0.98, 0.7),
        borderaxespad=0,
        frameon=True
    )
    fig.add_artist(legend1)

    # Size legend (Gene Ratio)
    example_ratios = [0.25, 0.6, 1.0]
    size_handles = [
        plt.scatter([], [], s=ratio * 500, color='grey', edgecolor='black', label=f"{ratio:.2f}")
        for ratio in example_ratios
    ]
    legend2 = fig.legend(
        handles=size_handles,
        title='Gene Ratio',
        loc='center right',
        bbox_to_anchor=(0.98, 0.3),
        borderaxespad=0,
        frameon=True
    )
    fig.add_artist(legend2)


def dotplot_enrichment(xd_smart:pl.DataFrame):
    # Step 1: Determine order based on median FDR
    term_order_df = (
        xd_smart
        .group_by("termDescription")
        .agg(pl.median("falseDiscoveryRate").alias("medianFDR"))
        .sort("medianFDR")  # This sorts ascending (smallest FDR first)
        .with_row_index(name="term_order")  # assign integer order
    )
    # Step 2: Join back to original data
    xd_smart = xd_smart.join(term_order_df.select(["termDescription", "term_order"]), on="termDescription", how="left")
    # Step 3: Sort by the new 'term_order' column - remove descending=True to keep smallest FDR at top
    xd_smart = xd_smart.sort("term_order")  # Removed descending=True

    pdf = xd_smart.to_pandas()
    # Step 1: Preprocess
    pdf, direction_colors = prepare_data_for_plotting(pdf)

    # Calculate figure height
    n_terms = pdf['termDescription'].nunique()
    height_per_term = 0.5      # inches per term
    max_height = 50         # cap so it doesn't grow unbounded
    fig_height = min(n_terms * height_per_term, max_height)
    
    # Create figure with proper spacing for legends
    fig, ax = plt.subplots(figsize=(16, fig_height))  # Increased width
    
    # Main scatter
    sc = plot_enrichment_scatter(ax, pdf)
    
    # Adjust layout to make space for legends - leave more space on the right
    fig.subplots_adjust(left=0.3, right=0.70)  # Reduced right margin to leave more space
    
    # Colorbar
    sm = plt.cm.ScalarMappable(
        cmap='coolwarm',
        norm=plt.Normalize(vmin=pdf['-log10FDR'].min(), vmax=pdf['-log10FDR'].max())
    )
    sm.set_array([])
    
    cbar = plt.colorbar(sm, orientation='horizontal', ax=ax, fraction=0.03, pad=0.15, shrink=0.8)
    cbar.set_label('-log10(FDR)')
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.xaxis.set_label_position('top')

    # Legends
    add_custom_legends(fig, ax, direction_colors)
    
    # Ensure the figure is properly laid out
    plt.tight_layout()
    
    plt.show()


def altair_dotplot_enrichment(xd_smart: pl.DataFrame) -> alt.Chart:
    """Return an Altair chart equivalent of :func:`dotplot_enrichment`.

    Parameters
    ----------
    xd_smart:
        Long-format Polars DataFrame with the same columns expected by
        :func:`dotplot_enrichment`.

    Returns
    -------
    alt.Chart
        Interactive Altair scatter plot.
    """

    term_order_df = (
        xd_smart
        .group_by("termDescription")
        .agg(pl.median("falseDiscoveryRate").alias("medianFDR"))
        .sort("medianFDR")
        .with_row_index(name="term_order")
    )
    xd_smart = xd_smart.join(
        term_order_df.select(["termDescription", "term_order"]),
        on="termDescription",
        how="left",
    )
    xd_smart = xd_smart.sort("term_order")

    pdf = xd_smart.to_pandas()
    pdf, direction_colors = prepare_data_for_plotting(pdf)

    chart_height = min(pdf["termDescription"].nunique() * 20, 1000)

    direction_scale = alt.Scale(
        domain=list(direction_colors.keys()),
        range=list(direction_colors.values()),
    )

    chart = (
        alt.Chart(pdf, height=chart_height, width=600)
        .mark_circle(strokeWidth=2)
        .encode(
            x=alt.X("contrast:N", title="Contrast"),
            y=alt.Y(
                "termLabel:N",
                sort=pdf["termLabel"].tolist(),
                title="Gene Set Terms",
            ),
            size=alt.Size(
                "geneRatio:Q",
                legend=alt.Legend(title="Gene Ratio"),
                scale=alt.Scale(range=[0, 500]),
            ),
            color=alt.Color(
                "-log10FDR:Q",
                scale=alt.Scale(scheme="redblue"),
                legend=alt.Legend(title="-log10(FDR)"),
            ),
            stroke=alt.Color(
                "direction:N",
                scale=direction_scale,
                legend=alt.Legend(title="Direction"),
            ),
            tooltip=[
                "termID",
                "termDescription",
                "contrast",
                "geneRatio",
                "falseDiscoveryRate",
                "direction",
            ],
        )
    )

    return chart
