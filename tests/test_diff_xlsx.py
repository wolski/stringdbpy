import polars as pl
import pytest

from string_gsea.models.gsea_models import RankList, RankListCollection
from string_gsea.ranks_from_dea_xlsx import DiffXLSX


class TestDiffXLSXIntegration:
    """Integration tests for the DiffXLSX class."""

    @pytest.fixture
    def zip_path(self, mouse_xlsx_zip):
        return mouse_xlsx_zip

    @pytest.fixture
    def diff_xlsx(self, zip_path):
        return DiffXLSX(zip_path)

    def test_initialization(self, diff_xlsx):
        """Test that the DiffXLSX class initializes correctly."""
        assert diff_xlsx.dea_df is not None
        assert isinstance(diff_xlsx.dea_df, pl.DataFrame)
        assert diff_xlsx.dea_df.height > 0
        expected_columns = ["IDcolumn", "contrast", "statistic", "modelName", "nrPeptides"]
        for col in expected_columns:
            assert col in diff_xlsx.dea_df.columns

    def test_rank_dict_returns_collection(self, diff_xlsx):
        """Test that rank_dict returns a RankListCollection for one analysis type."""
        collection = diff_xlsx.rank_dict(which="pep_2_no_imputed")
        assert isinstance(collection, RankListCollection)
        assert collection.analysis == "pep_2_no_imputed"
        assert len(collection) > 0

    def test_rank_dict_content(self, diff_xlsx):
        """Test the content of the rank_dict output."""
        collection = diff_xlsx.rank_dict(which="pep_1")
        for rl in collection:
            assert isinstance(rl, RankList)
            assert rl.n_genes > 0
            assert rl.contrast in collection

    def test_rank_dict_all_analysis_types(self, diff_xlsx):
        """Test that each valid analysis type works."""
        for which in ("pep_1", "pep_1_no_imputed", "pep_2", "pep_2_no_imputed"):
            collection = diff_xlsx.rank_dict(which=which)
            assert collection.analysis == which
            assert len(collection) > 0

    def test_rank_dict_invalid_analysis_raises(self, diff_xlsx):
        with pytest.raises(ValueError, match="which must be one of"):
            diff_xlsx.rank_dict(which="invalid")

    def test_rank_dict_contrasts(self, diff_xlsx):
        """Test that the rank_dict method correctly extracts contrasts."""
        original_contrasts = diff_xlsx.dea_df.select(pl.col("contrast")).unique().to_series().to_list()
        collection = diff_xlsx.rank_dict(which="pep_1")
        for contrast in original_contrasts:
            assert contrast in collection.contrasts
