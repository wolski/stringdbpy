import polars as pl
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap


def separate_pivot_longer(df: pl.DataFrame) -> pl.DataFrame:
    xd = (
        df
        # 1) split into list-columns, replacing the originals
    .with_columns([
        pl.col("proteinIDs")        .str.split(",").alias("proteinIDs"),
        pl.col("proteinLabels")     .str.split(",").alias("proteinLabels"),
        pl.col("proteinInputLabels").str.split(",").alias("proteinInputLabels"),
        pl.col("proteinInputValues").str.split(",").alias("proteinInputValues"),
        pl.col("proteinRanks").str.split(",").alias("proteinRanks"),
    ])
    # 2) explode all four at once (DF now has 161 872 rows)
    .explode([
        "proteinIDs",
        "proteinLabels",
        "proteinInputLabels",
        "proteinInputValues",
        "proteinRanks",
    ])
    )
        # 3) cast the exploded string→float
    xd = xd.with_columns(
            pl.col("proteinInputValues").cast(pl.Float64),
            pl.col("proteinRanks").cast(pl.Float64)
    )
    return xd

def summarize_terms(xd: pl.DataFrame) -> pl.DataFrame:
    # Compute mean input value per term
    means = (
        xd
        .group_by(["contrast", "termID"])
        .agg(pl.col("proteinInputValues").mean().alias("meanInputValues"))
    )
    xd = xd.join(means, on=["contrast", "termID"] )
    # Filter
    xd = xd.filter(
        (pl.col("falseDiscoveryRate") < 0.05) &
        (pl.col("genesMapped") > 10)
    )
    # Unique contrasts and categories
    #contrast_list = xd["contrast"].unique().to_list()
    #category_list = xd["category"].unique().to_list()
    #categ = [c for c in category_list if c != "Publications"]
    return xd


def make_network(xdf: pl.DataFrame) -> nx.Graph:
    # 1) dedupe the edge table
    edges_pl = (
        xdf
        .select([
            "termID",
            "proteinLabels",
            "falseDiscoveryRate",
            "proteinInputValues",   # plural here
            "meanInputValues",
            "direction"
        ])
        .unique()
    )

    # 2) term nodes need a placeholder for proteinInputValues (plural)
    terms_pl = (
        edges_pl
        .select(["termID",
                "falseDiscoveryRate",
                "meanInputValues",
                "direction"])
        .unique()
        .with_columns([
            pl.lit("term").alias("nodeType"),
            # now plural, matching edges_pl
            pl.lit(None).cast(pl.Float64).alias("proteinInputValues")
        ])
        .rename({"termID": "name"})
        .select([
            "name",
            "nodeType",
            "falseDiscoveryRate",
            "meanInputValues",
            "proteinInputValues",
            "direction"
        ])
    )

    # 3) protein nodes carry the real proteinInputValues
    proteins_pl = (
        edges_pl
        .select(["proteinLabels", "proteinInputValues"])
        .unique()
        .with_columns([
            pl.lit("protein").alias("nodeType"),
            pl.lit(None).cast(pl.Float64).alias("falseDiscoveryRate"),
            pl.lit(None).cast(pl.Float64).alias("meanInputValues"),
            pl.lit(None).cast(pl.Utf8).alias("direction")
        ])
        .rename({
            "proteinLabels": "name"
            # proteinInputValues stays the same
        })
        .select([
            "name",
            "nodeType",
            "falseDiscoveryRate",
            "meanInputValues",
            "proteinInputValues",
            "direction"
        ])
    )

    # 4) now they have identical schemas → concat
    nodes_pl = pl.concat([terms_pl, proteins_pl], how="vertical")

    # 5) build the edge‐list
    edge_list_pl = (
        edges_pl
        .rename({"termID": "from", "proteinLabels": "to"})
        .select(["from", "to"])
    )

    # 6) convert to small Python structures
    nodes = nodes_pl.to_dicts()
    edges = edge_list_pl.to_dicts()

    # 7) build the graph
    G = nx.Graph()
    for node in nodes:
        name = node.pop("name")
        node_type = node.pop("nodeType")
        # all other keys (falseDiscoveryRate, meanInputValues, proteinInputValues)
        G.add_node(name, nodeType=node_type, **node)
    for e in edges:
        G.add_edge(e["from"], e["to"])

    return G




def assign_node_sizes(G):
    """
    Mutates G.nodes[n]['size']:
      - term  = -log10(falseDiscoveryRate) * 4
      - protein = constant = 3
    """
    for n, d in G.nodes(data=True):
        if d["nodeType"] == "term":
            # protect against FDR = 0
            fdr = max(d["falseDiscoveryRate"], 1e-300)
            d["size"] = -np.log10(fdr) * 4
        else:
            d["size"] = 3
    return G



def assign_node_colors(G, n_colors: int = 100):
    """
    Assigns G.nodes[n]['color'] for *all* nodes using a diverging
    green→gray→red palette with a pivot at 0:
      - term nodes use meanInputValues
      - protein nodes use proteinInputValues
    Both get semi-transparency of 0.5.
    """
    # build the palette
    cmap = LinearSegmentedColormap.from_list("diverge", ["green", "gray", "red"], N=n_colors)
    half = n_colors // 2

    def map_val_to_rgba(val, vmin, vmax):
        """Map a single value to an index in [0..n_colors-1], pivoting at zero."""
        if val < 0:
            # negative half: [vmin..0] → [0..half-1]
            frac = (val - vmin) / (0 - vmin) if vmin < 0 else 0.0
            idx  = int(round(frac * (half - 1)))
        else:
            # positive half: [0..vmax] → [half..n_colors-1]
            frac = val / vmax if vmax > 0 else 0.0
            idx  = half + int(round(frac * (half - 1)))
        return cmap(idx)

    # Process term nodes
    term_nodes = [n for n,d in G.nodes(data=True) if d["nodeType"] == "term"]
    term_vals  = np.array([G.nodes[n]["meanInputValues"] for n in term_nodes], dtype=float)
    tmin, tmax = term_vals.min(), term_vals.max()
    # Since pivot is 0, ensure ranges straddle zero?  If all positive, vmin→0 maps to center.
    for n, val in zip(term_nodes, term_vals):
        rgba = map_val_to_rgba(val, tmin, tmax)
        # give terms 50% opacity as well (or set 1.0 if you prefer)
        G.nodes[n]["color"] = (rgba[0], rgba[1], rgba[2], 0.5)

    # Process protein nodes
    prot_nodes = [n for n,d in G.nodes(data=True) if d["nodeType"] == "protein"]
    prot_vals  = np.array([G.nodes[n]["proteinInputValues"] for n in prot_nodes], dtype=float)
    pmin, pmax = prot_vals.min(), prot_vals.max()
    for n, val in zip(prot_nodes, prot_vals):
        rgba = map_val_to_rgba(val, pmin, pmax)
        G.nodes[n]["color"] = (rgba[0], rgba[1], rgba[2], 0.5)

    return G

def make_network_with_colors(xdf: pl.DataFrame) -> nx.Graph:
    G = make_network(xdf)
    G = assign_node_sizes(G)
    G = assign_node_colors(G)
    return G


def plot_network_graph(G, title: str):
    # layout
    pos = nx.kamada_kawai_layout(G)

    # split nodes by type
    term_nodes = [n for n,d in G.nodes(data=True) if d["nodeType"] == "term"]
    prot_nodes = [n for n,d in G.nodes(data=True) if d["nodeType"] == "protein"]

    # gather attributes
    term_sizes  = [G.nodes[n]["size"]    * 5   for n in term_nodes]   # bump up term sizes by 3×
    term_colors = [G.nodes[n]["color"]       for n in term_nodes]
    prot_sizes  = [G.nodes[n]["size"]        for n in prot_nodes]
    prot_colors = [G.nodes[n]["color"]       for n in prot_nodes]

    # draw proteins first, unlabeled
    nx.draw_networkx_nodes(
        G, pos,
        nodelist   = prot_nodes,
        node_size  = prot_sizes,
        node_color = prot_colors,
        alpha      = 0.7
    )
    # draw term nodes next, with labels
    nx.draw_networkx_nodes(
        G, pos,
        nodelist   = term_nodes,
        node_size  = term_sizes,
        node_color = term_colors,
        #edgecolors = "red",
        linewidths = 0.5
    )
    nx.draw_networkx_labels(
        G, pos,
        labels = {n: n for n in term_nodes},
        font_size = 6,
        font_color = "black",
        alpha      = 0.7
    )
    # draw edges underneath
    nx.draw_networkx_edges(
        G, pos,
        edge_color = "lightgray",
        alpha      = 0.5,
        width      = 0.5
    )

    plt.title(title)
    plt.axis("off")
    plt.tight_layout()
    plt.show()
