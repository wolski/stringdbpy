import polars as pl


def filter_by_FDR(
    xd: pl.DataFrame, FDR_threshold: float = 0.05, genes_mapped_threshold: int = 10
) -> pl.DataFrame:
    return xd.filter(
        (pl.col("falseDiscoveryRate") < FDR_threshold)
        & (pl.col("genesMapped") > genes_mapped_threshold)
    )


def select_top_terms(df: pl.DataFrame, max_terms: int = 100) -> pl.DataFrame:
    """Keep only the top max_terms terms per (contrast, category) by FDR."""
    return (
        df.sort(["contrast", "category", "falseDiscoveryRate"])
        .with_columns(
            pl.col("termID").rank(method="dense").over(["contrast", "category"]).alias("_rank")
        )
        .filter(pl.col("_rank") <= max_terms)
        .drop("_rank")
    )


def add_gene_ratio(df: pl.DataFrame) -> pl.DataFrame:
    return df.with_columns(
        (pl.col("genesMapped") / pl.col("genesInSet")).alias("geneRatio")
    )
    return df


# def explode_protein_columns(df: pl.DataFrame) -> pl.DataFrame:
# Rename function to be more descriptive
def explode_protein_columns(df: pl.DataFrame) -> pl.DataFrame:
    """
    Separate protein columns that contain comma-separated values into individual rows.

    This function takes a DataFrame with protein-related columns containing comma-separated values
    and splits them into separate rows, with one protein per row. The affected columns are:
    - proteinIDs
    - proteinLabels
    - proteinInputLabels
    - proteinInputValues
    - proteinRanks

    Args:
        df (pl.DataFrame): Input DataFrame with comma-separated protein columns

    Returns:
        pl.DataFrame: DataFrame with protein data split into separate rows
    """

    # Step 1a: protect commas in proteinLabels that are followed by 1–2 digits and another comma
    df_protected = df.with_columns(
        [
            pl.col("proteinLabels")
            .str.replace_all(r",(\d{1,2},)", r"§COMMA§$1")
            .alias("proteinLabels")
        ]
    )

    # Step 1b: split all relevant columns (no chaining)
    df_split = df_protected.with_columns(
        [
            pl.col("proteinIDs").str.split(",").alias("proteinIDs"),
            pl.col("proteinLabels").str.split(",").alias("proteinLabels"),
            pl.col("proteinInputLabels").str.split(",").alias("proteinInputLabels"),
            pl.col("proteinInputValues").str.split(",").alias("proteinInputValues"),
            pl.col("proteinRanks").str.split(",").alias("proteinRanks"),
        ]
    )

    # Step 1c: restore protected commas inside the split proteinLabels lists (no chaining)
    df_split = df_split.with_columns(
        pl.col("proteinLabels")
        .list.eval(pl.element().str.replace_all("§COMMA§", ","))
        .alias("proteinLabels")
    )

    # 2) explode all four at once (DF now has 161 872 rows)
    df_exploded = df_split.explode(
        [
            "proteinIDs",
            "proteinLabels",
            "proteinInputLabels",
            "proteinInputValues",
            "proteinRanks",
        ]
    )

    # 3) cast the exploded string→float
    xd = df_exploded.with_columns(
        pl.col("proteinInputValues").cast(pl.Float64),
        pl.col("proteinRanks").cast(pl.Float64),
    )
    return xd


def summarize_terms(xd: pl.DataFrame) -> pl.DataFrame:
    # Compute mean input value per term
    means = xd.group_by(["contrast", "termID"]).agg(
        pl.col("proteinInputValues").mean().alias("meanInputValues")
    )
    xd = xd.join(means, on=["contrast", "termID"])
    # Filter
    xd = xd.filter((pl.col("falseDiscoveryRate") < 0.05) & (pl.col("genesMapped") > 10))
    return xd


