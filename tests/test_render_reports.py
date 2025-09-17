import pytest
from pathlib import Path
from string_gsea.scripts.render_reports import render_quarto_docs, execute_quarto_command

class TestRenderQuartoDocs:
    @pytest.fixture
    def mock_execute_quarto(self, mocker):
        return mocker.patch('string_gsea.scripts.render_reports.execute_quarto_command', return_value=True)
     
    
    @pytest.fixture
    def data_dir(self, tmp_path):
        """Create a temporary data directory with mock files."""
        data_dir = tmp_path / "test_data"
        data_dir.mkdir()
        # Create mock files that will be found by glob
        results_file = data_dir / "test_string_gsea_results_long.xlsx"
        links_file = data_dir / "test_links.txt"
        results_file.touch()
        links_file.touch()
        return data_dir
    
    def test_render_quarto_docs_calls_execute_with_correct_args(self, mock_execute_quarto, data_dir):
        """Test that render_quarto_docs calls execute_quarto_command with expected args."""
        # Run the function
        render_quarto_docs(data_dir)
        
        # Check that execute_quarto_command was called with the right arguments
        mock_execute_quarto.assert_called_once()
        exec_args = mock_execute_quarto.call_args.kwargs
        assert exec_args['docs_path'].name == "reports"
        assert exec_args['output_dir'] == data_dir / "reports"
        assert exec_args['xlsx_file'].name == "test_string_gsea_results_long.xlsx"
        assert exec_args['links_file'].name == "test_links.txt"


class TestExecuteQuartoCommand:
    @pytest.fixture
    def mock_subprocess_run(self, mocker):
        """Mock subprocess.run for success case."""
        mock_run = mocker.patch('subprocess.run')
        # Configure the mock to return a successful result
        mock_process = mocker.MagicMock()
        mock_process.returncode = 0
        mock_process.stdout = "Success output"
        mock_process.stderr = ""
        mock_run.return_value = mock_process
        return mock_run
    
    def test_execute_quarto_command_success(self, mock_subprocess_run):
        """Test execute_quarto_command when subprocess.run succeeds."""
        # Test parameters
        docs_path = Path('/test/docs/path')
        output_dir = Path('/test/output/dir')
        xlsx_file = Path('/test/data/results.xlsx')
        links_file = Path('/test/data/links.txt')
        
        # Call the function
        result = execute_quarto_command(docs_path, output_dir, xlsx_file, links_file)
        
        # Verify subprocess.run was called with the expected command
        mock_subprocess_run.assert_called_once()
        args, kwargs = mock_subprocess_run.call_args
        cmd = args[0]
        
        # Verify the command structure
        assert cmd[0] == "quarto"
        assert cmd[1] == "render"
        assert cmd[2] == str(docs_path)
        assert "--to" in cmd and "html" in cmd
        assert "--output-dir" in cmd and str(output_dir) in cmd
        assert any(f"data_file:{xlsx_file}" in arg for arg in cmd)
        assert any(f"links_file:{links_file}" in arg for arg in cmd)
        
        # Verify subprocess.run was called with check=True and other expected kwargs
        assert kwargs.get('check') is True
        assert kwargs.get('capture_output') is True
        assert kwargs.get('text') is True
        
        # Verify the function returned True for success
        assert result is True