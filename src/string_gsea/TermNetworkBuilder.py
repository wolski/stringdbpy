import polars as pl
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge

class TermNetworkBuilder:
    def __init__(self, xd: pl.DataFrame, category: str = "SMART"):
        self.xd_s_smart = xd.filter(pl.col("category") == category).select(["contrast", "termID", "proteinLabels"])

    def compute_node_sizes(self) -> dict[str,int]:
        node_sizes_df = (
            self.xd_s_smart
            .group_by("termID")
            .agg(pl.col("proteinLabels").n_unique().alias("size"))
        )
        node_sizes = {r["termID"]: r["size"] for r in node_sizes_df.to_dicts()}
        return node_sizes


    def build_shared_counts(self) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
        """
        From the filtered SMART table, compute:
        • node_sizes: dict termID → total # unique proteins
        • within_df: edges within the same contrast (w)
        • cross_df: edges between different contrasts (w)
        • all_df:      union of within_df + cross_df (w)
        """
        
        # 2) Self‐join to count shared proteins
        sp = (
            self.xd_s_smart
            .join(self.xd_s_smart, on="proteinLabels", how="inner", suffix="_b")
            .group_by(["termID", "contrast", "termID_b", "contrast_b"])
            .agg(pl.count().alias("nr_proteins"))
            .rename({"termID":"termID_a", "contrast":"contrast_a"})
        )

        # 3) Split into within-contrast and cross-contrast
        within_df = (
            sp
            .filter(pl.col("contrast_a") == pl.col("contrast_b"))
            #.filter(pl.col("termID_a") != pl.col("termID_b"))
            .group_by(["termID_a", "termID_b"])
            .agg(pl.col("nr_proteins").sum().alias("w"))
        )

        cross_df = (
            sp
            .filter(pl.col("contrast_a") != pl.col("contrast_b"))
            #.filter(pl.col("termID_a") != pl.col("termID_b"))
            .group_by(["termID_a", "termID_b"])
            .agg(pl.col("nr_proteins").sum().alias("w"))
        )

        
        all_df = (
            sp
            .group_by(["termID_a", "termID_b"])
            .agg(pl.col("nr_proteins").sum().alias("w"))
        )
        return within_df, cross_df, all_df


    def build_contrast_counts(self) -> tuple[dict[str,dict[str,int]], list[str]]:
        """
        Compute for each termID, and for each contrast, how many unique proteins.
        Returns:
        • contrast_counts: { termID: { contrast_name: count, ... }, ... }
        • contrasts:       sorted list of all contrast names
        """
        tc_counts_df = (
            self.xd_s_smart
            .group_by(["termID", "contrast"])
            .agg(pl.col("proteinLabels").n_unique().alias("count"))
        )
        # build nested dict
        contrast_counts: dict[str, dict[str, int]] = {}
        for row in tc_counts_df.to_dicts():
            contrast_counts.setdefault(row["termID"], {})[row["contrast"]] = row["count"]
        # unique sorted contrast names
        contrasts = sorted(tc_counts_df.select("contrast").unique().to_series().to_list())

        return contrast_counts, contrasts

