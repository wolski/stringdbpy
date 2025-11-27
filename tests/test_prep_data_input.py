import pytest
import zipfile
from string_gsea.scripts.render_reports import prepare_data_input


class TestPrepareDataInput:
    """Test the 4 core scenarios for prepare_data_input function."""

    @pytest.fixture
    def sample_dir(self, tmp_path):
        """Create a sample directory with test files."""
        data_dir = tmp_path / "sample_data"
        data_dir.mkdir()
        (data_dir / "test_string_gsea_results_long.xlsx").write_text("xlsx content")
        (data_dir / "test_links.txt").write_text("links content")
        return data_dir

    @pytest.fixture
    def sample_zip(self, tmp_path):
        """Create a sample zip file."""
        zip_path = tmp_path / "sample_data.zip"
        with zipfile.ZipFile(zip_path, "w") as zf:
            zf.writestr("test_string_gsea_results_long.xlsx", "xlsx content")
            zf.writestr("test_links.txt", "links content")
        return zip_path

    def test_zip_input(self, sample_zip, tmp_path):
        """Test: input zip with specified output."""
        output_dir = tmp_path / "output"
        result = prepare_data_input(sample_zip, output_dir)

        assert result == output_dir
        assert (output_dir / "test_string_gsea_results_long.xlsx").exists()
        assert (output_dir / "test_links.txt").exists()

    def test_dir_input(self, sample_dir, tmp_path):
        """Test: input dir with specified output."""
        output_dir = tmp_path / "output"
        result = prepare_data_input(sample_dir, output_dir)

        assert result == output_dir
        assert (output_dir / "test_string_gsea_results_long.xlsx").exists()
        assert (output_dir / "test_links.txt").exists()
