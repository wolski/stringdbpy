import pytest

from string_gsea.gsea_utilities import get_rank_files, write_rank_files
from string_gsea.models.gsea_models import RankList, RankListCollection


def test_get_rank_files_reads_from_DE_zip(dummy_session_rnk_zip):
    """Tests reading .rnk files from the DE_*.zip archive."""
    if not dummy_session_rnk_zip.exists():
        pytest.fail(f"Test zip file not found at {dummy_session_rnk_zip}")

    collection = get_rank_files(dummy_session_rnk_zip)

    assert isinstance(collection, RankListCollection)
    assert collection.analysis == "from_rnk"
    assert len(collection) > 0, f"No .rnk files found or read from {dummy_session_rnk_zip}."

    for rl in collection:
        assert isinstance(rl, RankList)
        assert rl.n_genes > 0, f"RankList for {rl.contrast}.rnk is empty."


def test_write_rank_files(tmp_path):
    """Tests the standalone write_rank_files function."""
    ranks = RankListCollection(
        analysis="pep_1",
        rank_lists=[
            RankList(contrast="contrast_A", entries={"G1": 1.5, "G2": -0.3}),
            RankList(contrast="contrast_B", entries={"G3": 2.1}),
        ],
    )
    outdir = tmp_path / "output"
    outdir.mkdir()

    result = write_rank_files(ranks, outdir)

    assert result == outdir
    rnk_files = list(outdir.glob("**/*.rnk"))
    assert len(rnk_files) == 2

    content = (outdir / "pep_1" / "contrast_A.rnk").read_text()
    assert "G1\t1.5" in content
    assert "G2\t-0.3" in content
