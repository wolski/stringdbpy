import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
from matplotlib import cm
import polars as pl
from string_gsea.TermNetworkBuilder import TermNetworkBuilder
from string_gsea.network import plot_network_graph_plotly
import plotly.graph_objs as go


class TermNetworkPlotter:
    def __init__(self,
                 node_sizes: dict[str,int],
                 contrast_counts: dict[str,dict[str,int]],
                 contrasts: list[str],
                 max_radius: float = 0.1):
        self.node_sizes      = node_sizes
        self.contrast_counts = contrast_counts
        self.contrasts       = contrasts
        self.max_radius      = max_radius
        self._max_size       = max(node_sizes.values())
         # pick a qualitative colormap with enough distinct colours
        n = len(contrasts)
        if n == 1:
            # Single contrast: just pick first color from tab10
            cmap = cm.get_cmap("tab10")
            self.colors = [cmap(0)]
        elif n <= 10:
            cmap = cm.get_cmap("tab10")
            self.colors = [cmap(i/(n-1)) for i in range(n)]
        elif n <= 20:
            cmap = cm.get_cmap("tab20")
            self.colors = [cmap(i/(n-1)) for i in range(n)]
        else:
            # Fallback to default property cycle if >20 contrasts
            palette = plt.rcParams['axes.prop_cycle'].by_key()['color']
            # Repeat palette to cover all contrasts
            palette *= ((n // len(palette)) + 1)
            self.colors = palette[:n]

    def _prepare_graph(
        self,
        edge_df: pl.DataFrame,
        thresh: int = 2,
        include_all: bool = False
    ) -> tuple[nx.Graph, pl.DataFrame, list[str]]:
        """
        Shared helper: filter edges at `thresh`, determine node set,
        and build NetworkX graph with those nodes and edges.

        If include_all=True, include all self.node_sizes keys as nodes;
        otherwise only nodes that appear in filtered edges.
        """
        # filter edges once
        edges_sub = edge_df.filter(pl.col("w") >= thresh)

        # decide which nodes to add
        if include_all:
            nodes = list(self.node_sizes.keys())
        else:
            nodes = (
                edges_sub
                .select(["termID_a","termID_b"])
                .melt(value_name="node")
                .select("node")
                .unique()
                .to_series()
                .to_list()
            )

        # build the graph of those nodes + edges
        G = nx.Graph()
        G.add_nodes_from(nodes)
        for a,b,w in edges_sub.select(["termID_a","termID_b","w"]).iter_rows():
            if a != b:
                G.add_edge(a, b, weight=w)
        return G, edges_sub, nodes
    

    def compute_full_layout(
        self,
        edge_df: pl.DataFrame,
        thresh: int = 2,
        k: float = 0.8,
        seed: int = 0
    ) -> dict[str,tuple[float,float]]:
        """
        Compute and store a spring layout over the *full* term set.
        """
        G, _, _ = self._prepare_graph(edge_df, thresh=thresh, include_all=True)
        # compute layout once (includes isolates)
        self.fixed_pos = nx.spring_layout(G, k=k, seed=seed)
        return self.fixed_pos

    

    def _draw_pies(self,
                   ax: plt.Axes,
                   pos: dict[str,tuple[float,float]],
                   node_sizes: dict[str,int],
                   max_size: int) -> None:
        """Draw pie charts at each node position."""
        for term,(x,y) in pos.items():
            counts = [self.contrast_counts.get(term,{}).get(c,0)
                      for c in self.contrasts]
            total = sum(counts)
            if total == 0:
                continue
            r = self.max_radius * (node_sizes[term]/max_size)
            start = 0.0
            for frac,color in zip(counts,self.colors):
                angle = 360*frac/total
                wedge = Wedge(
                    center=(x,y),
                    r=r,
                    theta1=start,
                    theta2=start+angle,
                    facecolor=color,
                    edgecolor="white",
                    linewidth=0.5,
                    zorder=2
                )
                ax.add_patch(wedge)
                start += angle

    def draw_panel(self,
                   ax: plt.Axes,
                   edge_df,
                   thresh: int = 2,
                   title: str = "",
                   use_fixed_layout: bool = False) -> None:
        """
        Draw a network panel with nodes sized/pie colored and edges filtered by thresh.
        """
        # use helper to get graph, filtered edges, and node list
        G, edges_sub, nodes = self._prepare_graph(
            edge_df, thresh=thresh, include_all=False
        )
        
        # If no nodes, just set title and return
        if not nodes:
            ax.set_title(title, fontsize=10, fontweight='bold')
            ax.axis('off')
            ax.text(0.5, 0.5, 'No terms meet threshold', 
                   ha='center', va='center', transform=ax.transAxes,
                   fontsize=10, color='gray')
            return
        
        # 3) build temp sizes dict for those nodes
        temp_sizes = {n: self.node_sizes[n] for n in nodes}
        temp_max   = max(temp_sizes.values()) if temp_sizes else 1

                 
        # 5) layout
        if use_fixed_layout and self.fixed_pos is not None:
            full_pos = self.fixed_pos
            pos = {n: full_pos[n] for n in nodes}
            xs, ys = zip(*(pos[n] for n in nodes))
            pad = 0.1 * max(max(xs)-min(xs), max(ys)-min(ys))
            ax.set_xlim(min(xs)-pad, max(xs)+pad)
            ax.set_ylim(min(ys)-pad, max(ys)+pad)
        else:
            pos = nx.spring_layout(G, k=0.8, seed=0)

        # 6) draw edges
        if G.edges():
            ec = nx.draw_networkx_edges(
                G, pos,
                width=[edge_attrs['weight']/5 for _, _, edge_attrs in G.edges(data=True)],
                edge_color="#888888",
                ax=ax
            )
            if ec is not None:  # Check if ec is not None before setting zorder
                ec.set_zorder(1)

        # 7) draw pies
        self._draw_pies(ax, pos, temp_sizes, temp_max)

        # 8) labels
        labs = nx.draw_networkx_labels(G, pos, font_size=6, ax=ax)
        for txt in labs.values():
            txt.set_zorder(3)

        ax.set_title(title)
        ax.axis("off")

    
    @staticmethod
    def get_figure_legend(one_contrast:bool = True, category:str = "SMART", thresh:int = 1) -> str:
        if one_contrast:
            text = f"""Network representations of {category} co-enrichment across contrasts (T ≥ {thresh} proteins).
                    Node size is proportional to the total number of proteins annotated to each term; edge width reflects the number of shared proteins.
                    """
        else:
            text = (
                    f"""Network representations of {category} co-enrichment across contrasts (T ≥ {thresh} proteins).\n
                    (A) Full network: all term–term pairs that share at least {thresh} proteins either within a single contrast or across distinct contrasts.
                    Node size is proportional to the total number of proteins annotated to each term; edge width reflects the number of shared proteins.\n
                    (B) Cross‐contrast subnetwork: links between two categories (terms) where one term comes from contrast X and the other comes from contrast Y, but they share the same proteins,
                    highlighting pathway relationships that persist between contrasts.\n
                    (C) Within‐contrast subnetwork: links between two terms (categories) that come from the same contrast.
                    revealing context-specific functional links that may dissolve upon perturbation.\n"""
            )
        return text

    def draw_legend_panel(
        self,
        ax: plt.Axes,
        title: str = "Contrast",
        loc: str = "center"
    ) -> None:
        """
        Draw a standalone legend of contrasts and their colors in the given Axes.
        """
        import matplotlib.patches as mpatches
        handles = [
            mpatches.Patch(facecolor=col, edgecolor="black", label=ctr)
            for ctr, col in zip(self.contrasts, self.colors)
        ]
        ax.legend(
            handles=handles,
            title=title,
            loc=loc,
            frameon=False,
            fontsize=8,
            title_fontsize=9
        )
        ax.axis("off")
    


def plot_network(xd, category:str="SMART", contrast:str=None, thresh:int=3, use_fixed_layout:bool=True):
    if contrast is None:
        nb = TermNetworkBuilder(xd, category=category)
    else:
        nb = TermNetworkBuilder(xd.filter(pl.col("contrast") == contrast), category=category)
    
    within_df, cross_df, all_df = nb.build_shared_counts()
    contrast_counts, contrasts = nb.build_contrast_counts()
    node_sizes = nb.compute_node_sizes()

    plotter = TermNetworkPlotter(
        node_sizes      = node_sizes,
        contrast_counts = contrast_counts,
        contrasts       = contrasts,
        max_radius      = 0.1
    )
    plotter.compute_full_layout(all_df, thresh=1)

    # Set higher DPI for better resolution
    plt.rcParams['figure.dpi'] = 300  # Default is usually 100

    if contrast is None:
        fig, axes = plt.subplots(4, 1, figsize=(7,18))
        plotter.draw_panel(axes[0], all_df,   thresh=thresh, title=f"A) Full (T≥{thresh})",
                            use_fixed_layout=use_fixed_layout)
        plotter.draw_panel(axes[2], cross_df, thresh=thresh, title=f"C) Cross (T≥{thresh})",
                            use_fixed_layout=use_fixed_layout)
        plotter.draw_panel(axes[1], within_df,thresh=thresh, title=f"D) Within (T≥{thresh})",
                            use_fixed_layout=use_fixed_layout)
        plotter.draw_legend_panel(axes[3], title="Contrast", loc="center")
    else:
        fig, axes = plt.subplots(1, 1, figsize=(7,7))
        plotter.draw_panel(axes, all_df, thresh=thresh, title=f"Network for {contrasts[0]}, threshold = {thresh}",
                          use_fixed_layout=use_fixed_layout)
    
    plt.tight_layout()
    plt.show()


