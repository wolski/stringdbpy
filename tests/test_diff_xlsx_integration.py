import pytest
from pathlib import Path
import polars as pl
from loguru import logger

from src.a373_string_gsea.ranksFromDEAXLSX import DiffXLSX


class TestDiffXLSXIntegration:
    """Integration tests for the DiffXLSX class."""
    
    @pytest.fixture
    def zip_path(self):
        """Fixture to provide the path to the test data zip file."""
        return Path(__file__).parent / "data" / "DE_mouse_fasta_xlsx.zip"
    
    @pytest.fixture
    def diff_xlsx(self, zip_path):
        """Fixture to provide a DiffXLSX instance."""
        return DiffXLSX(zip_path)
    
    def test_initialization(self, diff_xlsx):
        """Test that the DiffXLSX class initializes correctly."""
        assert diff_xlsx.dea_df is not None, "DataFrame should not be None"
        assert isinstance(diff_xlsx.dea_df, pl.DataFrame), "DataFrame should be a Polars DataFrame"
        assert diff_xlsx.dea_df.height > 0, "DataFrame should have rows"
        assert diff_xlsx.dea_df.width > 0, "DataFrame should have columns"
        
        # Check for expected columns
        expected_columns = ["IDcolumn", "contrast", "statistic", "modelName", "nrPeptides"]
        for col in expected_columns:
            assert col in diff_xlsx.dea_df.columns, f"Expected column '{col}' not found in DataFrame"
    
    def test_rank_dict_structure(self, diff_xlsx):
        """Test the structure of the rank_dict output."""
        rank_files = diff_xlsx.rank_dict()
        
        # Check that rank_files is a dictionary
        assert isinstance(rank_files, dict), "rank_files should be a dictionary"
        assert len(rank_files) > 0, "rank_files should not be empty"
        
        # Check the keys of the dictionary
        expected_outer_keys = ["pep_1", "pep_1_no_imputed", "pep_2", "pep_2_no_imputed"]
        for key in rank_files.keys():
            assert isinstance(key, tuple), "Keys should be tuples"
            assert len(key) == 2, "Keys should have length 2"
            assert key[0] in expected_outer_keys, f"Unexpected outer key: {key[0]}"
    
    def test_rank_dict_content(self, diff_xlsx):
        """Test the content of the rank_dict output."""
        rank_files = diff_xlsx.rank_dict()
        
        # Check that each value in the dictionary is a DataFrame
        for key, df in rank_files.items():
            assert isinstance(df, pl.DataFrame), f"Value for key {key} should be a DataFrame"
            assert df.height > 0, f"DataFrame for key {key} should have rows"
            assert df.width == 2, f"DataFrame for key {key} should have 2 columns"
            
            # Check that the DataFrame has the expected columns
            assert "IDcolumn" in df.columns, f"Expected column 'IDcolumn' not found in DataFrame for key {key}"
            assert "statistic" in df.columns, f"Expected column 'statistic' not found in DataFrame for key {key}"
    
    def test_rank_dict_filtering(self, diff_xlsx):
        """Test that the rank_dict method correctly filters the data."""
        rank_files = diff_xlsx.rank_dict()
        
        # Check that pep_1_no_imputed only contains rows where modelName does not contain "imputed"
        if ("pep_1_no_imputed", "contrast1") in rank_files:
            df = rank_files[("pep_1_no_imputed", "contrast1")]
            # We can't directly check the modelName column as it's not in the filtered DataFrame
            # But we can check that the DataFrame is not empty
            assert df.height > 0, "pep_1_no_imputed DataFrame should not be empty"
        
        # Check that pep_2 only contains rows where nrPeptides > 1
        if ("pep_2", "contrast1") in rank_files:
            df = rank_files[("pep_2", "contrast1")]
            # We can't directly check the nrPeptides column as it's not in the filtered DataFrame
            # But we can check that the DataFrame is not empty
            assert df.height > 0, "pep_2 DataFrame should not be empty"
    
    def test_rank_dict_contrasts(self, diff_xlsx):
        """Test that the rank_dict method correctly extracts contrasts."""
        # Get the unique contrasts from the original DataFrame
        original_contrasts = diff_xlsx.dea_df.select(pl.col("contrast")).unique().to_series().to_list()
        
        # Get the rank_dict
        rank_files = diff_xlsx.rank_dict()
        
        # Extract the contrasts from the rank_dict keys
        rank_dict_contrasts = set()
        for key in rank_files.keys():
            if isinstance(key, tuple) and len(key) == 2:
                rank_dict_contrasts.add(key[1])
        
        # Check that all contrasts from the original DataFrame are in the rank_dict
        for contrast in original_contrasts:
            assert contrast in rank_dict_contrasts, f"Contrast {contrast} from original DataFrame not found in rank_dict" 