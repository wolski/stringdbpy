import pytest
import polars as pl
import networkx as nx
from pathlib import Path
from matplotlib.figure import Figure

from string_gsea import network


TEST_DATA_DIR = Path(__file__).parent / 'data'

@pytest.fixture
def long_results_df():
    # Load a small sample of the Excel file for testing
    file_path = TEST_DATA_DIR/"dummy_out/WU_abcd_GSEA/from_rnk/WUabcd_string_gsea_results_long.xlsx"
    # Only read the first 100 rows for speed
    df = pl.read_excel(file_path)
    cont = "Bait_NCP_pUbT12T14_results.tsv"
    cat = "Monarch"
    df = df.filter((pl.col("contrast")==cont) & (pl.col("category")==cat))
    return df

def test_separate_pivot_longer(long_results_df):
    xd = network.separate_pivot_longer(long_results_df)
    assert isinstance(xd, pl.DataFrame)
    # Check that columns are exploded to lists
    assert "proteinIDs" in xd.columns
    assert xd.height > long_results_df.height

def test_summarize_terms(long_results_df):
    xd = network.separate_pivot_longer(long_results_df)
    summarized = network.summarize_terms(xd)
    assert isinstance(summarized, pl.DataFrame)
    # Should have filtered for FDR and genesMapped
    assert all(summarized["falseDiscoveryRate"] < 0.05)
    assert all(summarized["genesMapped"] > 10)

def test_make_network_with_colors(long_results_df):
    xd = network.separate_pivot_longer(long_results_df)
    summarized = network.summarize_terms(xd)
    G = network.make_network_with_colors(summarized)
    assert isinstance(G, nx.Graph)
    # Check that nodes have color and size attributes
    for n, d in G.nodes(data=True):
        assert "color" in d
        assert "size" in d

def test_plot_network_graph_runs(long_results_df):
    xd = network.separate_pivot_longer(long_results_df)
    summarized = network.summarize_terms(xd)
    G = network.make_network_with_colors(summarized)
    # Just check that plotting runs without error
    plt = network.plot_network_graph(G, title="Test Plot")
    assert plt is not None
    print(type(plt))
    # check type of plt
    assert isinstance(plt, Figure)

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
    xd = network.separate_pivot_longer(long_results_df)
    summarized = network.summarize_terms(xd)
    G = network.make_network_with_colors(summarized)
    widget = network.interactive_cytoscape(G)
    from ipycytoscape import CytoscapeWidget
    assert isinstance(widget, CytoscapeWidget)

def test_plot_network_graph_plotly(long_results_df):
    xd = network.separate_pivot_longer(long_results_df)
    summarized = network.summarize_terms(xd)
    G = network.make_network_with_colors(summarized)
    fig = network.plot_network_graph_plotly(G, title="Test Plotly")
    import plotly.graph_objs as go
    assert isinstance(fig, go.Figure)

    