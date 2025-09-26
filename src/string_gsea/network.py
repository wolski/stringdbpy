import polars as pl
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from ipycytoscape import CytoscapeWidget
import networkx as nx
import plotly.graph_objs as go
import matplotlib.pyplot as plt
import networkx as nx

def filter_by_FDR(xd: pl.DataFrame,
 FDR_threshold: float = 0.05,
 genes_mapped_threshold: int = 10) -> pl.DataFrame:
    return xd.filter(
        (pl.col("falseDiscoveryRate") < FDR_threshold) &
        (pl.col("genesMapped") > genes_mapped_threshold)
    )


def add_gene_ratio(df: pl.DataFrame) -> pl.DataFrame:
    return df.with_columns(
        (pl.col("genesMapped") / pl.col("genesInSet")).alias("geneRatio")
    )
    return df


#def explode_protein_columns(df: pl.DataFrame) -> pl.DataFrame:
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
    df_protected = df.with_columns([
        pl.col("proteinLabels").str.replace_all(r",(\d{1,2},)", r"§COMMA§$1").alias("proteinLabels")
    ])

    # Step 1b: split all relevant columns (no chaining)
    df_split = df_protected.with_columns([
        pl.col("proteinIDs").str.split(",").alias("proteinIDs"),
        pl.col("proteinLabels").str.split(",").alias("proteinLabels"),
        pl.col("proteinInputLabels").str.split(",").alias("proteinInputLabels"),
        pl.col("proteinInputValues").str.split(",").alias("proteinInputValues"),
        pl.col("proteinRanks").str.split(",").alias("proteinRanks"),
    ])

    # Step 1c: restore protected commas inside the split proteinLabels lists (no chaining)
    df_split = df_split.with_columns(
        pl.col("proteinLabels").list.eval(
            pl.element().str.replace_all("§COMMA§", ",")
        ).alias("proteinLabels")
    )

    if False:
            # Step 2: Check if all split lists have equal length
        df_with_lengths = df_split.with_columns([
            pl.col("proteinIDs").list.len().alias("proteinIDs_len"),
            pl.col("proteinLabels").list.len().alias("proteinLabels_len"),
            pl.col("proteinInputLabels").list.len().alias("proteinInputLabels_len"),
            pl.col("proteinInputValues").list.len().alias("proteinInputValues_len"),
            pl.col("proteinRanks").list.len().alias("proteinRanks_len"),
        ])

        # Find rows where lengths don't match
        problematic_rows = df_with_lengths.filter(
            (pl.col("proteinIDs_len") != pl.col("proteinLabels_len")) |
            (pl.col("proteinIDs_len") != pl.col("proteinInputLabels_len")) |
            (pl.col("proteinIDs_len") != pl.col("proteinInputValues_len")) |
            (pl.col("proteinIDs_len") != pl.col("proteinRanks_len"))
        )


    # 2) explode all four at once (DF now has 161 872 rows)
    df_exploded = df_split.explode([
        "proteinIDs",
        "proteinLabels",
        "proteinInputLabels",
        "proteinInputValues",
        "proteinRanks",
    ])
    
        # 3) cast the exploded string→float
    xd = df_exploded.with_columns(
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
        description = data.get("termDescription", "")
        return (
            f"<b>Term:</b> {node_id}<br/>"
            f"<b>Direction:</b> {direction}<br/>"
            f"<b>Mean:</b> {mean_val:.3f}<br/>"
            f"<b>FDR:</b> {fdr:.3f}<br/>"
            f"<b>Description:</b> {description}"
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

def bipartite_hybrid_layout(
    G,
    terms: list,
    proteins: list,
    left_x: float = -1.0,
    right_x: float = +1.0,
    input_attr: str = "proteinInputValues",
    desc: bool = True
):
    """
    1) Sort proteins by G.nodes[p][input_attr] descending (or ascending if desc=False).
    2) Lay terms in input-sorted order? No—keep original term order for now.
    3) Compute term-barycenters over that protein order, re-sort terms.
    4) Compute protein-barycenters over the new term order, re-sort proteins.
    5) Lay both sides evenly in y.
    """
    # (A) fetch values & initial sort of proteins
    prot_vals = {p: G.nodes[p].get(input_attr, 0.0) for p in proteins}
    prots_init = sorted(proteins, key=lambda p: prot_vals[p], reverse=desc)

    # (B) initial term positions (in their given order)
    #    we'll postpone final term sort until after step C
    T = len(terms)
    yL0 = np.linspace(1.0, -1.0, T)
    term_pos = {t: (left_x, yL0[i]) for i,t in enumerate(terms)}

    # (C) compute term-barycenters over the initial protein order
    #     use the y from an evenly spaced prot_init layout
    P0 = len(prots_init)
    yP0 = np.linspace(1.0, -1.0, P0)
    prot_pos0 = {p: (right_x, yP0[i]) for i,p in enumerate(prots_init)}

    term_bary = {}
    for t in terms:
        nbr = [p for p in G[t] if p in prot_pos0]
        if nbr:
            term_bary[t] = np.mean([prot_pos0[p][1] for p in nbr])
        else:
            term_bary[t] = 0.0

    sorted_terms = sorted(terms, key=lambda t: term_bary[t], reverse=True)

    # (D) re-lay out terms in new barycenter order
    yL1 = np.linspace(1.0, -1.0, len(sorted_terms))
    term_pos = {t: (left_x, yL1[i]) for i,t in enumerate(sorted_terms)}

    # (E) now compute protein-barycenters over the new term positions
    prot_bary = {}
    for p in proteins:
        nbr = [t for t in G[p] if t in term_pos]
        if nbr:
            prot_bary[p] = np.mean([term_pos[t][1] for t in nbr])
        else:
            prot_bary[p] = 0.0

    # (F) final protein order: start from the input-sorted list, 
    #     then sub-sort by barycenter to refine crossings
    sorted_prots = sorted(
        prots_init,
        key=lambda p: prot_bary[p],
        reverse=True
    )

    # (G) lay out proteins in that final order
    yP1 = np.linspace(1.0, -1.0, len(sorted_prots))
    prot_pos = {p: (right_x, yP1[i]) for i,p in enumerate(sorted_prots)}

    # (H) merge & return
    return {**term_pos, **prot_pos}


def bipartite_barycenter_layout(
    G,
    terms,      # list of "term" node IDs (original set)
    prots,      # list of "protein" node IDs (original set)
    left_x=-1.0,
    right_x=+1.0
):
    # 1) Start by laying out terms in their given order
    L = len(terms)
    yL = np.linspace( 1, -1, L )
    term_pos = { t: (left_x, yL[i]) for i,t in enumerate(terms) }

    # 2) First barycenter pass: sort proteins by avg-y of their term neighbours
    prot_bary = {}
    for p in prots:
        nbr = [t for t in G[p] if t in term_pos]
        prot_bary[p] = np.mean([term_pos[t][1] for t in nbr]) if nbr else 0.0
    # descending so highest-bary (towards top) appear at top
    sorted_prots = sorted(prots, key=lambda p: prot_bary[p], reverse=True)

    # 3) Lay out proteins in that new order
    P = len(sorted_prots)
    yP = np.linspace( 1, -1, P )
    prot_pos = { p: (right_x, yP[i]) for i,p in enumerate(sorted_prots) }

    # 4) Second pass: re-sort terms by barycenters of these protein positions
    term_bary = {}
    for t in terms:
        nbr = [p for p in G[t] if p in prot_pos]
        term_bary[t] = np.mean([prot_pos[p][1] for p in nbr]) if nbr else 0.0
    sorted_terms = sorted(terms, key=lambda t: term_bary[t], reverse=True)

    # 5) Re-lay out terms in that final order
    yL2 = np.linspace( 1, -1, len(sorted_terms) )
    term_pos = { t: (left_x, yL2[i]) for i,t in enumerate(sorted_terms) }

    # 6) Return merged dict
    return {**term_pos, **prot_pos}


def plot_network_graph_plotly(G, title: str, layout: str = "bipartite"):
    # 1) layout
        # 1) pick layout
    if layout == "kamada_kawai":
        pos = nx.kamada_kawai_layout(G)
    elif layout == "circular":
        pos = nx.circular_layout(G)
    elif layout == "bipartite":
        terms = [n for n,d in G.nodes(data=True) if d["nodeType"]=="term"]
        prots = [n for n,d in G.nodes(data=True) if d["nodeType"]=="protein"]
        pos   = bipartite_hybrid_layout(G, terms, prots)
    else:
        raise ValueError(f"Unknown layout: {layout!r}")
    
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
        build_tooltip(n, G.nodes[n]).replace("<br/>", "<br>")
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
    term_hover = [
        build_tooltip(n, G.nodes[n]).replace("<br/>", "<br>")
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
