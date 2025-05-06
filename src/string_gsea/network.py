import polars as pl
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from ipycytoscape import CytoscapeWidget
import networkx as nx
import plotly.graph_objs as go


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
            "direction",
            "termDescription"  # Added termDescription
        ])
        .unique()
    )

    # 2) term nodes need a placeholder for proteinInputValues (plural)
    terms_pl = (
        edges_pl
        .select(["termID",
                "falseDiscoveryRate",
                "meanInputValues",
                "direction",
                "termDescription"])
        .unique()
        .with_columns([
            pl.lit("term").alias("nodeType"),
            pl.lit(None).cast(pl.Float64).alias("proteinInputValues")
        ])
        .rename({"termID": "name"})
        .select([
            "name",
            "nodeType",
            "falseDiscoveryRate",
            "meanInputValues",
            "proteinInputValues",
            "direction",
            "termDescription"
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
            pl.lit(None).cast(pl.Utf8).alias("direction"),
            pl.lit(None).cast(pl.Utf8).alias("termDescription")
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
            "direction",
            "termDescription"
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
    graph = nx.Graph()
    for node in nodes:
        name = node.pop("name")
        node_type = node.pop("nodeType")
        # all other keys (falseDiscoveryRate, meanInputValues, proteinInputValues)
        graph.add_node(name, nodeType=node_type, **node)
    for e in edges:
        graph.add_edge(e["from"], e["to"])

    return graph


def assign_node_sizes(G):
    """
    Mutates G.nodes[n]['size']:
      - term  = -log10(falseDiscoveryRate) * 4
      - protein = constant = 3
    """
    for n, d in G.nodes(data=True):
        if d["nodeType"] == "term":
            # protect against FDR = 0
            fdr = max(d["falseDiscoveryRate"], 1e-5)
            d["size"] = -np.log10(fdr) * 3
        else:
            d["size"] = 3
    return G


def assign_node_colors(G, n_colors: int = 100):
    """
    Assigns for each node in G:
      - G.nodes[n]['color']       → RGBA fill (diverging green→gray→red, pivot=0)
      - G.nodes[n]['borderColor'] → 'yellow' / 'blue' / 'green' / 'gray' (by direction)
    """
    # 1) build the diverging fill palette
    cmap = LinearSegmentedColormap.from_list("diverge", ["green", "gray", "red"], N=n_colors)
    half = n_colors // 2

    def map_val_to_rgba(val, vmin, vmax):
        """Map a value to an index in [0..n_colors-1], pivoting at zero, and return the RGBA tuple."""
        if val < 0:
            # negative half: map [vmin..0] → [0..half-1]
            frac = (val - vmin) / (0 - vmin) if (vmin < 0) else 0.0
            idx  = int(round(frac * (half - 1)))
        else:
            # positive half: map [0..vmax] → [half..n_colors-1]
            frac = val / vmax if (vmax > 0) else 0.0
            idx  = half + int(round(frac * (half - 1)))
        return cmap(idx)

    # 2) map directions → border-colors
    direction_color_map = {
        'top':       'yellow',
        'bottom':    'blue',
        'both_ends': 'green'
    }

    # 3) Term nodes: compute fill & border
    term_nodes = [n for n,d in G.nodes(data=True) if d["nodeType"] == "term"]
    term_vals  = np.array([G.nodes[n]["meanInputValues"] for n in term_nodes], dtype=float)
    tmin, tmax = term_vals.min(), term_vals.max()

    for n, val in zip(term_nodes, term_vals):
        rgba = map_val_to_rgba(val, tmin, tmax)
        G.nodes[n]["color"] = (rgba[0], rgba[1], rgba[2], 0.5)

        # direction → borderColor (with safe default)
        direction = G.nodes[n].get("direction", "").replace(" ", "_").lower()
        G.nodes[n]["borderColor"] = direction_color_map.get(direction, "gray")

    # 4) Protein nodes: compute fill & assign default borderColor
    prot_nodes = [n for n,d in G.nodes(data=True) if d["nodeType"] == "protein"]
    prot_vals  = np.array([G.nodes[n]["proteinInputValues"] for n in prot_nodes], dtype=float)
    pmin, pmax = prot_vals.min(), prot_vals.max()

    for n, val in zip(prot_nodes, prot_vals):
        rgba = map_val_to_rgba(val, pmin, pmax)
        G.nodes[n]["color"]       = (rgba[0], rgba[1], rgba[2], 0.5)
        G.nodes[n]["borderColor"] = "#cccccc"  # or whatever default

    return G


def make_network_with_colors(xdf: pl.DataFrame) -> nx.Graph:
    G = make_network(xdf)
    G = assign_node_sizes(G)
    G = assign_node_colors(G)
    return G


import matplotlib.pyplot as plt
import networkx as nx

def plot_network_graph(G, title: str):
    # 1) compute layout
    pos = nx.kamada_kawai_layout(G)

    # 2) split nodes by type
    term_nodes = [n for n,d in G.nodes(data=True) if d["nodeType"] == "term"]
    prot_nodes = [n for n,d in G.nodes(data=True) if d["nodeType"] == "protein"]

    # 3) gather attributes
    term_sizes        = [G.nodes[n]["size"] * 5 for n in term_nodes]   # scale up terms
    term_face_colors  = [G.nodes[n]["color"]   for n in term_nodes]
    term_edge_colors  = [G.nodes[n].get("borderColor", "gray") for n in term_nodes]

    prot_sizes        = [G.nodes[n]["size"]    for n in prot_nodes]
    prot_face_colors  = [G.nodes[n]["color"]   for n in prot_nodes]

    # 4) draw edges first so they sit underneath
    nx.draw_networkx_edges(
        G, pos,
        edge_color = "lightgray",
        alpha      = 0.5,
        width      = 0.5
    )

    # 5) draw protein nodes (unlabeled, semi‐transparent)
    nx.draw_networkx_nodes(
        G, pos,
        nodelist    = prot_nodes,
        node_size   = prot_sizes,
        node_color  = prot_face_colors,
        alpha       = 0.7,
        edgecolors  = None  # or a default if you prefer borders on proteins
    )

    # 6) draw term nodes with colored borders
    nx.draw_networkx_nodes(
        G, pos,
        nodelist    = term_nodes,
        node_size   = term_sizes,
        node_color  = term_face_colors,
        edgecolors  = term_edge_colors,
        linewidths  = 1.0,    # thick enough to see the colored border
        alpha       = 0.5
    )

    # 7) add labels for terms
    nx.draw_networkx_labels(
        G, pos,
        labels     = {n: n for n in term_nodes},
        font_size  = 6,
        font_color = "black",
        alpha      = 0.8
    )

    # 8) finish styling
    plt.title(title)
    plt.axis("off")
    plt.tight_layout()
    fig = plt.gcf()  # Get the current figure
    return fig



def _rgba_to_css(rgba):
    r, g, b, a = rgba
    return f"rgba({int(r*255)},{int(g*255)},{int(b*255)},{a:.2f})"


def build_tooltip(node_id: str, data: dict) -> str:
    """
    Generate the HTML tooltip for a Cytoscape node.

    Args:
      node_id: the node’s identifier (usually your label for term nodes).
      data:    the node’s attribute dict, must contain:
               - "nodeType" (either "term" or "protein")
               - for terms: "direction" and "meanInputValues"
               - for proteins: "proteinInputValues"

    Returns:
      A HTML string suitable for Cytoscape’s tooltip-text.
    """
    if data["nodeType"] == "term":
        direction = data.get("direction", "unknown")
        mean_val  = data.get("meanInputValues", float("nan"))
        fdr       = data.get("falseDiscoveryRate", float("nan"))
        return (
            f"<b>Term:</b> {node_id}<br/>"
            f"<b>Direction:</b> {direction}<br/>"
            f"<b>Mean:</b> {mean_val:.3f}<br/>"
            f"<b>FDR:</b> {fdr:.3f}"
        )
    else:
        prot_val = data.get("proteinInputValues", float("nan"))
        return (
            f"<b>Protein:</b> {node_id}<br/>"
            f"<b>Value:</b> {prot_val:.3f}"
        )


def interactive_cytoscape(G, layout="cose", width="100%", height="600px"):
    nodes = []
    for n, d in G.nodes(data=True):
        is_term = (d["nodeType"] == "term")
        size    = d["size"] * (2 if is_term else 1)
        color   = _rgba_to_css(d["color"])
        label   = n if is_term else ""
        tip     = build_tooltip(n, d)  # your existing tooltip logic

        # now just read borderColor out of the node data
        border_color = d.get("borderColor", "black")

        nodes.append({
            "data": {
                "id":          n,
                "label":       label,
                "tooltip":     tip,
                "nodeType":    d["nodeType"],
                "size":        size,
                "color":       color,
                "borderColor": border_color
            },
            # you can keep classes for other styling if you like
            "classes": "term" if is_term else "protein"
        })

    edges = [{"data": {"source": u, "target": v}} for u, v in G.edges()]

    cw = CytoscapeWidget()
    cw.layout.width  = width
    cw.layout.height = height
    cw.graph.add_graph_from_json({"nodes": nodes, "edges": edges})
    cw.set_layout(name=layout)

    cw.set_style([
        # 1) Base node rule (unchanged)
        {
            "selector": "node",
            "style": {
                "width":            "data(size)",
                "height":           "data(size)",
                "background-color": "data(color)",
                "content":          "data(label)",
                "font-size":        "8px",
                "text-valign":      "center",
                "text-halign":      "center",
                "overlay-padding":  "6px",
                "tooltip-text":     "data(tooltip)",
                "tooltip-show":     "yes"
            }
        },
        # 2) Protein nodes (unchanged)
        {
            "selector": "node.protein",
            "style": {
                "content":            "",
                "background-opacity": 0.5
            }
        },
        # 3) Term nodes: use data(borderColor)
        {
            "selector": "node.term",
            "style": {
                "background-color":   "data(color)",
                "background-opacity": 1,
                "border-width":       "2px",
                "border-color":       "data(borderColor)"
            }
        },
        # 4) Edges (unchanged)
        {
            "selector": "edge",
            "style": {
                "line-color":   "#cccccc",
                "line-opacity": 0.5,
                "width":        1
            }
        }
    ])

    return cw



def plot_network_graph_plotly(G, title: str):
    # 1) layout
    pos = nx.kamada_kawai_layout(G)

    def to_css_rgba(rgba):
        r, g, b, a = rgba
        return f"rgba({int(r*255)},{int(g*255)},{int(b*255)},{a})"

    # 2) make the edge trace
    edge_x, edge_y = [], []
    for u, v in G.edges():
        x0, y0 = pos[u]; x1, y1 = pos[v]
        edge_x += [x0, x1, None]; edge_y += [y0, y1, None]
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        mode="lines",
        line=dict(width=0.5, color="lightgray"),
        hoverinfo="none"
    )

    # 3) proteins
    prot_nodes  = [n for n,d in G.nodes(data=True) if d["nodeType"]=="protein"]
    prot_x      = [pos[n][0] for n in prot_nodes]
    prot_y      = [pos[n][1] for n in prot_nodes]
    prot_sizes  = [G.nodes[n]["size"]    for n in prot_nodes]
    prot_colors = [to_css_rgba(G.nodes[n]["color"]) for n in prot_nodes]
    # use plain <br> here
    prot_hover  = [
        f"Protein: {n}<br>Value: {G.nodes[n]['proteinInputValues']:.3f}"
        for n in prot_nodes
    ]
    prot_trace = go.Scatter(
        x=prot_x, y=prot_y,
        mode="markers",
        marker=dict(
            size=prot_sizes,
            color=prot_colors,
            opacity=0.7,
            line=dict(width=0)    # no border
        ),
        hoverinfo="text",
        hovertext=prot_hover,
        # force a neutral white box with black text
        hoverlabel=dict(
            bgcolor="white",
            bordercolor="black",
            font=dict(color="black")
        ),
        showlegend=False
    )

    # 4) terms
    term_nodes       = [n for n,d in G.nodes(data=True) if d["nodeType"]=="term"]
    term_x           = [pos[n][0] for n in term_nodes]
    term_y           = [pos[n][1] for n in term_nodes]
    term_sizes       = [G.nodes[n]["size"]*5 for n in term_nodes]
    term_face_colors = [to_css_rgba(G.nodes[n]["color"]) for n in term_nodes]
    term_edge_colors = [G.nodes[n].get("borderColor","gray")        for n in term_nodes]
    term_hover       = [
        # again use <br> for line breaks
        f"Term: {n}<br>"
        f"Direction: {G.nodes[n].get('direction','unknown')}<br>"
        f"Mean: {G.nodes[n]['meanInputValues']:.3f}"
        for n in term_nodes
    ]
    term_trace = go.Scatter(
        x=term_x, y=term_y,
        mode="markers+text",
        marker=dict(
            size=term_sizes,
            color=term_face_colors,
            opacity=0.5,
            line=dict(width=1, color=term_edge_colors)
        ),
        text=term_nodes,
        textposition="top center",
        textfont=dict(size=8, color="black"),
        hoverinfo="text",
        hovertext=term_hover,
        hoverlabel=dict(
            bgcolor="white",
            bordercolor="black",
            font=dict(color="black")
        ),
        showlegend=False
    )

    # 5) assemble
    fig = go.Figure(data=[edge_trace, prot_trace, term_trace])
    fig.update_layout(
        title=title,
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        plot_bgcolor="white",
        margin=dict(l=20, r=20, t=40, b=20)
    )
    return fig
