
import polars as pl
import seaborn as sns
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage


def pivot_to_wide(df: pl.DataFrame) -> pl.DataFrame:
    """
    1. Pivot the dataframe to wide format
    2. Convert values to binary (0/1) based on presence
    3. Group by contrast and category
    """

    df_sorted = df.sort("proteinInputValues", descending=False)
    # First pivot to wide format
    wide_df = df_sorted.pivot(
        values="proteinInputValues",
        index="proteinLabels",
        on="termID"
    )
    
    return wide_df

def make_nested_dict(df: pl.DataFrame) -> dict:
    nested = {}
    
    # get all contrasts
    contrasts = df["contrast"].unique().to_list()
    for c in contrasts:
        # slice down to this contrast
        df_c = df.filter(pl.col("contrast") == c)
        
        # get all categories in that contrast
        categories = df_c["category"].unique().to_list()
        cat_dict = {}
        for cat in categories:
            # slice down to this category
            df_cat = df_c.filter(pl.col("category") == cat)
            df_cat = df_cat.drop(["contrast", "category"])
            cat_dict[cat] = pivot_to_wide(df_cat)
        nested[c] = cat_dict
    
    return nested

def convert_to_binary(df: pl.DataFrame, index_col: str = "proteinLabels", to_boolean: bool = False) -> pl.DataFrame:
    dtype = pl.Boolean if to_boolean else pl.Int8
    binary_df = df.with_columns([
        pl.col(c)
        .is_not_null()
        .cast(dtype)
        .alias(c)
        for c in df.columns
        if c != index_col
    ])
    return binary_df


def plot_term_distance_heatmap(fc: pl.DataFrame) -> None:
    
    
    binary_df = convert_to_binary(fc)
    binary_pdf = (
        binary_df.to_pandas()
        .set_index("proteinLabels")
    )
    fc_pdf = (
        fc.to_pandas()
        .set_index("proteinLabels")
    )

    # 1) Compute linkages from the binary matrix
    #    (rows = terms, cols = proteins)
    row_dist = pdist(binary_pdf.values, metric="jaccard")
    col_dist = pdist(binary_pdf.values.T, metric="jaccard")

    row_link = linkage(row_dist, method="average")
    col_link = linkage(col_dist, method="average")

    # 2) Call clustermap on the *fold‐change* matrix,
    #    but re‐use the precomputed linkages for ordering
    g = sns.clustermap(
        fc_pdf,
        row_linkage=row_link,
        col_linkage=col_link,
        cmap="RdBu_r",          # e.g. a divergent cmap for fold‐changes
        center=0,             # or 0.0, depending on what “no‐change” is
        figsize=(12, 12),
        xticklabels=True,      # still hide the crowded protein names
        yticklabels=False
    )
    return g
