import upsetplot
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from scipy.stats import gaussian_kde
import pandas as pd
from upsetplot import from_indicators, UpSet

def make_upset(df_binary: pl.DataFrame, df_values: pl.DataFrame, max_category = 25,max_subset_rank = 35) -> None:
    df_values = df_values.select(pl.col("proteinLabels"), pl.col("proteinInputValues")).unique()
    df = df_values.join(df_binary, on="proteinLabels", how="right")
    df_binary = df.to_pandas()

    # Get domain columns (all columns except proteinLabels and proteinInputValues)
    non_domain_cols = ['proteinLabels', 'proteinInputValues']
    domain_cols = [col for col in df_binary.columns if col not in non_domain_cols]

    # compute total hits per domain, pick top N
    counts = df_binary[domain_cols].sum(axis=0).sort_values(ascending=False)
    topN   = counts.head(max_category).index.tolist()

    # rebuild only those columns
    df_small = df_binary[topN + non_domain_cols]
    #pdfx     = df_small.set_index(topN)
    df_small = df_small[df_small[topN].any(axis=1)]
    pdfx = df_small.set_index(list(topN))
    # limit the number of elements in the intersection plot
    upset = upsetplot.UpSet(pdfx, subset_size="count",intersection_plot_elements=3,max_subset_rank= max_subset_rank)

    upset.add_catplot(value="proteinInputValues", kind="strip", color="blue")
    return upset


def plot_single_ridge(
    ax,
    x_vals: np.ndarray,
    gene_ratio: float,
    norm_fdr: float,
    direction: str,
    term_id: str,
    description: str,
    x_min: float,
    x_max: float,
    desc_max_chars: int = 60,      # NEW: max display length
    desc_fontsize: float = 6.0     # NEW: smaller font
):
    """
    Draw a single ridge plot for a gene set enrichment term.

    Args:
        ax: Matplotlib axis to draw on
        x_vals (np.ndarray): Input values for the proteins in this term
        gene_ratio (float): Ratio of genes in list vs term size, between 0-1
        norm_fdr (float): Normalized false discovery rate, between 0-1
        direction (str): Direction of enrichment, one of 'top', 'bottom', or 'both ends'
        term_id (str): ID of the enrichment term
        description (str): Description of the enrichment term
        x_min (float): Minimum x value for plotting range
        x_max (float): Maximum x value for plotting range
    """

    """
    Draw one ridge on `ax` with:
      - KDE of x_vals from x_min→x_max
      - fill color = gene_ratio (0→white, 1→red)
      - border color = {'top':'yellow','bottom':'blue','both ends':'green'}[direction]
      - a tiny bar at the right whose height = norm_fdr (0→1)
      - y=0 baseline, term_id label, and centered description
    """
    # compute sweep
    xi = np.linspace(x_min, x_max, 200)
    yi = gaussian_kde(x_vals)(xi)
    
    # colors
    fill_col   = cm.Greens(gene_ratio)
    border_col = {"top":"yellow", "bottom":"blue", "both ends":"green"}[direction]
    
    # draw
    ax.fill_between(xi, yi, color=fill_col, alpha=0.8)
    ax.plot(xi, yi, color=border_col, lw=1.5)
    ax.axhline(0, color="k", lw=0.5)
    
    # little FDR bar
    x_range   = x_max - x_min
    bar_off   = 0.05 * x_range
    pad       = 0.02 * x_range
    bar_x     = x_max + bar_off/2
    bar_w     = bar_off
    ax.bar(bar_x, norm_fdr, width=bar_w, bottom=0, color="grey", alpha=0.7)
    
    # fix x-limits so both data & bar are in view
    ax.set_xlim(x_min - pad, x_max + bar_off + pad)
    
    
    # annotate
    ax.set_yticks([])
    ax.set_ylabel(term_id, rotation=0, labelpad=50, va="center", fontsize=9)

    if len(description) > desc_max_chars:
        display_desc = description[: desc_max_chars - 3] + "..."
    else:
        display_desc = description

    ax.text(0.5, 0.5, display_desc,
            transform=ax.transAxes,
            ha="center", va="center",
            fontsize=desc_fontsize)

def plot_term_ridges(df_long, ridge_height=0.4, ridge_width=6, left_margin=0.2, hspace=0.4):
    # — prepare pandas df —
    pdf = (
        df_long
        .select([
            "termID","termDescription",
            "proteinInputValues","geneRatio",
            "falseDiscoveryRate","direction"   # renamed from 'border'
        ])
        .filter(pl.col("proteinInputValues").is_not_null())
        .to_pandas()
    )
    
    # order by median
    term_order = (
        pdf.groupby("termID")["proteinInputValues"]
           .median()
           .sort_values(ascending=False)
           .index
           .tolist()
    )
    cmap = cm.get_cmap("Reds")
    # global min/max
    x_min, x_max = pdf["proteinInputValues"].min(), pdf["proteinInputValues"].max()
    # precompute normalized FDR per term
    neglog = -np.log10(pdf["falseDiscoveryRate"])
    fdr_max = neglog.max()
    
    # build canvas
    n = len(term_order)
    fig, axes = plt.subplots(n, 1,
                             figsize=(ridge_width, ridge_height*n),
                             sharex=False, dpi=300)
    if n == 1:
        axes = [axes]
    
    # loop
    for ax, tid in zip(axes, term_order):
        sub = pdf[pdf.termID == tid]
        plot_single_ridge(
            ax               = ax,
            x_vals           = sub["proteinInputValues"].values,
            gene_ratio       = sub["geneRatio"].iat[0],
            norm_fdr         = (-np.log10(sub["falseDiscoveryRate"].iat[0]) / fdr_max),
            direction        = sub["direction"].iat[0],
            term_id          = tid,
            description      = sub["termDescription"].iat[0],
            x_min            = x_min,
            x_max            = x_max,
        )
    
    # only bottom axis shows x‐ticks
    for ax in axes[:-1]:
        ax.tick_params(labelbottom=False)
    axes[-1].set_xlabel("proteinInputValues")
    
    plt.subplots_adjust(hspace=hspace, left=left_margin, bottom=0.1)
    sns.despine(left=True, bottom=True)
    return fig


def make_upset_contrasts_terms(xd:pl.DataFrame, category:str="SMART"):
    """
    Make an upset plot of the contrasts and terms

    Args:
        xd: The dataframe to plot
        category: The category of the terms to plot
    """
    xd_s_smart = (
        xd
        .filter(pl.col("category") == category)
        .select(["termID", "contrast"])
        .unique()
    )

    # 2) Convert to a pandas indicator matrix: index=termID, columns=contrast, True if present
    df_pd = xd_s_smart.to_pandas()
    indicator = pd.crosstab(df_pd['termID'], df_pd['contrast']).astype(bool)

    # 3) Build the UpSet data structure
    upset_data = from_indicators(indicator.columns.tolist(), indicator)

    # 4) Plot
    plt.figure(figsize=(8,5))
    u = UpSet(
        upset_data,
        show_counts=True,
        sort_by='degree',        # orders intersections by size
        sort_categories_by=None,
    )

    u.plot()
    plt.suptitle(f"{category} term Presence Across Contrasts")
    plt.show()
