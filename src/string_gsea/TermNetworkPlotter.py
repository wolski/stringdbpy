import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
from matplotlib import cm
import polars as pl

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


    def compute_full_layout(self, edge_df, thresh: int = 2, k: float = 0.8, seed: int = 0):
        """
        Build the graph from edge_df at threshold and compute a fixed spring layout.
        Stores positions in self.fixed_pos.
        """
        G = nx.Graph()
        for term, sz in self.node_sizes.items():
            G.add_node(term, size=sz)
        for termA, termB, w in edge_df.select(["termID_a","termID_b","w"]).iter_rows():
            if w >= thresh:
                G.add_edge(termA, termB, weight=w)
        G.remove_nodes_from(list(nx.isolates(G)))
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
        # 1) filter edges
        edges_sub = edge_df.filter(pl.col("w") >= thresh)

        # 2) participating nodes via Polars
        nodes = (
            edges_sub
            .select(["termID_a", "termID_b"])
            .melt(value_name="node")
            .select("node")
            .unique()
            .to_series()
            .to_list()
        )

        # 3) build temp sizes dict for those nodes
        temp_sizes = {n: self.node_sizes[n] for n in nodes}
        temp_max   = max(temp_sizes.values()) if temp_sizes else 1

        # 4) build graph
        G = nx.Graph()
        G.add_nodes_from(nodes)
        for a, b, w in edges_sub.select(["termID_a","termID_b","w"]).iter_rows():
            if a != b:
                G.add_edge(a, b, weight=w)
            
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
        ec = nx.draw_networkx_edges(
            G, pos,
            width=[edge_attrs['weight']/5 for _, _, edge_attrs in G.edges(data=True)],
            edge_color="#888888",
            ax=ax
        )
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
    def get_figure_legend():
        text = """
            Figure:. Network representations of SMART‐term co-enrichment across two contrasts (T ≥ 2 proteins).
            (A) Full network: all term–term pairs that share at least two proteins either within a single contrast or across both experiments. Node size is proportional to the total number of proteins annotated to each term; edge width reflects the number of shared proteins.
            (B) Cross‐contrast subnetwork: only those edges supported by the same proteins in contrast 1 and contrast 2, highlighting pathway relationships that persist between conditions.
            (C) Within‐contrast subnetwork: only edges observed within a single experiment, revealing context-specific functional links that may dissolve upon perturbation.
            """
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

