import polars as pl
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
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

def plot_enrichment_scatter(ax, pdf):
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

    ax.margins(x=0.2)
    ax.set_xlabel('Contrast')
    ax.set_ylabel('Gene Set Terms')
    ax.set_xticks(range(len(pdf['contrast'].unique())))
    ax.set_xticklabels(pdf['contrast'].unique(), rotation=45, ha='right')


def add_custom_legends(fig, ax, direction_colors):
    # Border legend
    border_legend = [
        mlines.Line2D([], [], marker='o', linestyle='None', markersize=10,
                      markerfacecolor='white', markeredgecolor=color, label=label)
        for label, color in direction_colors.items()
    ]

    legend1 = ax.legend(
        handles=border_legend,
        title="Direction",
        loc='lower left',
        bbox_to_anchor=(-0.55, -0.20),
        frameon=True
    )

    # Size legend
    example_ratios = [0.25, 0.6, 1.0]
    size_scaler = lambda x: x * 500
    size_legend = [
        plt.scatter([], [], s=size_scaler(r), color='grey', edgecolor='black', label=f"{r:.2f}")
        for r in example_ratios
    ]

    legend2 = ax.legend(
        handles=size_legend,
        title="Gene Ratio",
        loc='lower left',
        bbox_to_anchor=(-0.55, -0.40),
        frameon=True
    )

    fig.add_artist(legend1)
    fig.add_artist(legend2)

def dotplot_enrichment(xd_smart:pl.DataFrame):
    # Step 1: Determine order based on median FDR
    term_order_df = (
        xd_smart
        .group_by("termDescription")
        .agg(pl.median("falseDiscoveryRate").alias("medianFDR"))
        .sort("medianFDR")
        .with_row_index(name="term_order")  # assign integer order
    )
    # Step 2: Join back to original data
    xd_smart = xd_smart.join(term_order_df.select(["termDescription", "term_order"]), on="termDescription", how="left")
    # Step 3: Sort by the new 'term_order' column
    xd_smart = xd_smart.sort("term_order", descending=True)

    pdf = xd_smart.to_pandas()
    # Step 1: Preprocess
    pdf, direction_colors = prepare_data_for_plotting(pdf)

    # Step 2: Create figure
    fig, ax = plt.subplots(figsize=(10, len(pdf['termDescription'].unique()) * 0.5))

    # Step 3: Main scatterplot
    plot_enrichment_scatter(ax, pdf)

    # Step 4: Add legends
    add_custom_legends(fig, ax, direction_colors)

    # Step 5: Add colorbar
    sm = plt.cm.ScalarMappable(
        cmap='coolwarm',
        norm=plt.Normalize(vmin=pdf['-log10FDR'].min(), vmax=pdf['-log10FDR'].max())
    )
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('-log10(FDR)')

    plt.tight_layout()
    plt.show()