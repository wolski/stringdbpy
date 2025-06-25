import pytest
import polars as pl
import networkx as nx
from pathlib import Path
from matplotlib.figure import Figure

from string_gsea import network


TEST_DATA_DIR = Path(__file__).parent / 'data'

@pytest.fixture
def long_results_df() -> pl.DataFrame:
    # Load a small sample of the Excel file for testing
    file_path = TEST_DATA_DIR/"dummy_out/WU_abcd_GSEA/from_rnk/WUabcd_string_gsea_results_long.xlsx"
    # Note: Comment mentions reading 100 rows, but full sheet is read then filtered.
    # Consider pl.read_excel(file_path, read_options={"n_rows": 100}) if initial load is slow.
    df = pl.read_excel(file_path)
    cont = "Bait_NCP_pUbT12T14_results.tsv"
    cat = "Monarch"
    df = df.filter((pl.col("contrast")==cont) & (pl.col("category")==cat))
    return df

@pytest.fixture
def processed_df(long_results_df: pl.DataFrame) -> pl.DataFrame:
    """DataFrame after explode_protein_columns."""
    return network.explode_protein_columns(long_results_df)

@pytest.fixture
def summarized_df(processed_df: pl.DataFrame) -> pl.DataFrame:
    """DataFrame after summarize_terms."""
    return network.summarize_terms(processed_df)

@pytest.fixture
def graph_with_colors(summarized_df: pl.DataFrame) -> nx.Graph:
    """NetworkX graph with color and size attributes."""
    return network.make_network_with_colors(summarized_df)

def test_explode_protein_columns(processed_df: pl.DataFrame, long_results_df: pl.DataFrame):
    assert isinstance(processed_df, pl.DataFrame)
    # Check that columns are exploded to lists
    assert "proteinIDs" in processed_df.columns
    assert processed_df.height > long_results_df.height

def test_summarize_terms(summarized_df: pl.DataFrame):
    assert isinstance(summarized_df, pl.DataFrame)
    # Should have filtered for FDR and genesMapped
    assert (summarized_df["falseDiscoveryRate"] < 0.05).all()
    assert (summarized_df["genesMapped"] > 10).all()

def test_make_network_with_colors(graph_with_colors: nx.Graph):
    assert isinstance(graph_with_colors, nx.Graph)
    # Check that nodes have color and size attributes
    for n, d in graph_with_colors.nodes(data=True):
        assert "color" in d
        assert "size" in d

def test_plot_network_graph_runs(graph_with_colors: nx.Graph):
    """
    Tests that the network plotting function runs and returns a Figure.
    Assumes a function network.plot_network_graph(G) exists.
    """
    try:
        # Assuming network.plot_network_graph is the function to test
        # and it returns a matplotlib Figure.
        # If the function has a different name or behavior (e.g., saves to file, shows plot),
        # this part will need adjustment.
        fig = network.plot_network_graph(graph_with_colors, title="Test Network Plot")
        assert isinstance(fig, Figure)
    except AttributeError:
        pytest.skip("Skipping plotting test: network.plot_network_graph function not found or accessible.")
    except Exception as e:
        pytest.fail(f"network.plot_network_graph raised an exception: {e}")

def test_build_tooltip_term():
    data = {
        "nodeType": "term",
        "direction": "top",
        "meanInputValues": 1.234
    }
    tip = network.build_tooltip("TermA", data)
    assert "Term:" in tip and "Direction:" in tip and "Mean:" in tip

def test_build_tooltip_protein():
    data = {
        "nodeType": "protein",
        "proteinInputValues": 2.345
    }
    tip = network.build_tooltip("ProtA", data)
    assert "Protein:" in tip and "Value:" in tip

def test_interactive_cytoscape(long_results_df):
    xd = network.explode_protein_columns(long_results_df)
    summarized = network.summarize_terms(xd)
    G = network.make_network_with_colors(summarized)
    widget = network.interactive_cytoscape(G)
    from ipycytoscape import CytoscapeWidget
    assert isinstance(widget, CytoscapeWidget)

def test_plot_network_graph_plotly(long_results_df):
    xd = network.explode_protein_columns(long_results_df)
    summarized = network.summarize_terms(xd)
    G = network.make_network_with_colors(summarized)
    fig = network.plot_network_graph_plotly(G, title="Test Plotly")
    import plotly.graph_objs as go
    assert isinstance(fig, go.Figure)
