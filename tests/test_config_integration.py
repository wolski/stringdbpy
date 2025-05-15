import os
import shutil
import tempfile
from pathlib import Path
import pytest
from unittest.mock import patch, MagicMock

from string_gsea.config import get_configuration, write_initial_configuration


@pytest.fixture
def temp_config_dir():
    """Create a temporary directory for testing configuration."""
    # Create a temporary directory
    temp_dir = tempfile.mkdtemp()
    
    # Mock the home directory path
    with patch('string_gsea.config.Path.home') as mock_home:
        mock_home.return_value = Path(temp_dir)
        
        # Mock the APPDATA environment variable for Windows tests
        with patch.dict('os.environ', {'APPDATA': temp_dir}):
            # Define the expected config directory and file
            config_dir = Path(temp_dir) / '.config' / 'string_gsea'
            config_file = config_dir / 'config.toml'
            
            yield {
                'temp_dir': temp_dir,
                'config_dir': config_dir,
                'config_file': config_file
            }
    
    # Clean up after the test
    shutil.rmtree(temp_dir)


@pytest.fixture
def mock_api_response():
    """Mock API response for STRING-DB API key request."""
    return [{"api_key": "test_api_key_123", "note": "This key will be activated within 30 minutes."}]


def test_write_initial_configuration(temp_config_dir, mock_api_response):
    """Test writing initial configuration with API key fetching."""
    # Mock the API response
    with patch('requests.get') as mock_get:
        mock_response = MagicMock()
        mock_response.json.return_value = mock_api_response
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response
        
        # Call the function with custom parameters
        fdr = 0.1
        caller_identity = "test.caller.com"
        
        config_path = write_initial_configuration(fdr=fdr, caller_identity=caller_identity)
        
        # Check that the config directory was created
        assert temp_config_dir['config_dir'].exists()
        
        # Check that the config file was created
        assert config_path.exists()
        assert config_path == temp_config_dir['config_file']
        
        # Verify the API was called
        mock_get.assert_called_once_with("https://version-12-0.string-db.org/api/json/get_api_key")
        
        # Read the config file and verify its contents
        with open(config_path, 'rb') as f:
            import tomli
            config = tomli.load(f)
        
        # Check that the config contains the expected values
        assert config['api_key'] == "test_api_key_123"
        assert config['fdr'] == fdr
        assert config['caller_identity'] == caller_identity
        assert config['ge_enrichment_rank_direction'] == 1


def test_get_configuration(temp_config_dir):
    """Test getting configuration from a file."""
    # First, create a config file
    temp_config_dir['config_dir'].mkdir(parents=True, exist_ok=True)
    
    # Create a sample config file
    sample_config = {
        'api_key': 'test_api_key_456',
        'fdr': 0.15,
        'ge_enrichment_rank_direction': -1,
        'caller_identity': 'test.caller.com'
    }
    
    with open(temp_config_dir['config_file'], 'wb') as f:
        import tomli_w
        tomli_w.dump(sample_config, f)
    
    # Get the configuration
    config = get_configuration()
    
    # Verify the configuration
    assert config['api_key'] == 'test_api_key_456'
    assert config['fdr'] == 0.15
    assert config['ge_enrichment_rank_direction'] == -1
    assert config['caller_identity'] == 'test.caller.com'


def test_get_configuration_missing_file(temp_config_dir):
    """Test getting configuration when the file doesn't exist."""
    # Ensure the config file doesn't exist
    if temp_config_dir['config_file'].exists():
        temp_config_dir['config_file'].unlink()
    
    # Expect a FileNotFoundError
    with pytest.raises(FileNotFoundError):
        get_configuration()


def test_get_configuration_missing_keys(temp_config_dir):
    """Test getting configuration when required keys are missing."""
    # Create a config file with missing keys
    temp_config_dir['config_dir'].mkdir(parents=True, exist_ok=True)
    
    # Create an incomplete config file
    incomplete_config = {
        'api_key': 'test_api_key_789',
        # Missing 'fdr' and 'caller_identity'
        'ge_enrichment_rank_direction': -1,
    }
    
    with open(temp_config_dir['config_file'], 'wb') as f:
        import tomli_w
        tomli_w.dump(incomplete_config, f)
    
    # Expect a ValueError about missing keys
    with pytest.raises(ValueError) as excinfo:
        get_configuration()
    
    # Check that the error message mentions the missing keys
    error_msg = str(excinfo.value)
    assert "missing required keys" in error_msg
    assert "fdr" in error_msg
    assert "caller_identity" in error_msg 