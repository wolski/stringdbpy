import os
import tomli
import tomli_w
import requests
import json
from pathlib import Path
from loguru import logger
from typing import Dict, Any, Optional

def get_configuration() -> Dict[str, Any]:
    """
    Get the configuration from a TOML file in a platform-independent way.
    
    The configuration file is located in:
    - Windows: %APPDATA%/string_gsea/config.toml
    - macOS/Linux: ~/.config/string_gsea/config.toml
    
    Returns:
        Dict[str, Any]: The configuration as a dictionary.
        
    Raises:
        FileNotFoundError: If the configuration file doesn't exist.
        ValueError: If the configuration file is invalid or missing required keys.
    """
    # Determine the config directory based on the operating system
    if os.name == 'nt':  # Windows
        config_dir = Path(os.environ.get('APPDATA', '')) / 'string_gsea'
    else:  # macOS and Linux
        config_dir = Path.home() / '.config' / 'string_gsea'
    
    # Path to the config file
    config_path = config_dir / 'config.toml'
    
    # If the config file doesn't exist, raise an exception
    if not config_path.exists():
        raise FileNotFoundError(
            f"Configuration file not found at {config_path}. "
            f"Please create a configuration file using write_initial_configuration() "
            f"and update it with your settings."
        )
    
    try:
        # Read the configuration file
        with open(config_path, 'rb') as f:
            config = tomli.load(f)
        
        # Check for required keys
        required_keys = ['api_key', 'fdr', 'ge_enrichment_rank_direction', 'caller_identity']
        missing_keys = [key for key in required_keys if key not in config]
        
        if missing_keys:
            raise ValueError(
                f"Configuration file is missing required keys: {', '.join(missing_keys)}. "
                f"Please update your configuration file at {config_path}."
            )
        
        return config
    except tomli.TOMLDecodeError as e:
        raise ValueError(f"Invalid TOML format in configuration file {config_path}: {e}")
    except Exception as e:
        raise ValueError(f"Error loading configuration from {config_path}: {e}")

def write_initial_configuration(fdr: float = 0.25, caller_identity: str = "www.fgcz.ch") -> Path:
    """
    Write an initial configuration file with the provided parameters and fetch an API key from STRING-DB.
    
    The configuration file is located in:
    - Windows: %APPDATA%/string_gsea/config.toml
    - macOS/Linux: ~/.config/string_gsea/config.toml
    
    Args:
        fdr: The false discovery rate threshold for GSEA analysis (default: 0.25)
        caller_identity: The caller identity for STRING-DB API (default: "www.fgcz.ch")
        
    Returns:
        Path: The path to the created configuration file.
        
    Raises:
        requests.RequestException: If there's an error fetching the API key from STRING-DB.
        ValueError: If the API key response is invalid.
    """
    # Determine the config directory based on the operating system
    if os.name == 'nt':  # Windows
        config_dir = Path(os.environ.get('APPDATA', '')) / 'string_gsea'
    else:  # macOS and Linux
        config_dir = Path.home() / '.config' / 'string_gsea'
    
    # Create the directory if it doesn't exist
    config_dir.mkdir(parents=True, exist_ok=True)
    
    # Path to the config file
    config_path = config_dir / 'config.toml'
    
    # Fetch API key from STRING-DB
    try:
        logger.info("Fetching API key from STRING-DB...")
        response = requests.get("https://version-12-0.string-db.org/api/json/get_api_key")
        response.raise_for_status()  # Raise an exception for HTTP errors
        
        api_key_data = response.json()
        if not api_key_data or not isinstance(api_key_data, list) or len(api_key_data) == 0:
            raise ValueError("Invalid response format from STRING-DB API")
            
        api_key = api_key_data[0].get("api_key")
        if not api_key:
            raise ValueError("No API key found in the response")
            
        logger.info(f"Successfully obtained API key: {api_key}")
        logger.info(f"Note: {api_key_data[0].get('note', 'No note provided')}")
    except requests.RequestException as e:
        logger.error(f"Error fetching API key from STRING-DB: {e}")
        raise
    except (json.JSONDecodeError, ValueError) as e:
        logger.error(f"Error parsing API key response: {e}")
        raise ValueError(f"Failed to parse API key response: {e}")
    
    # Configuration with fetched API key
    config = {
        'api_key': api_key,
        'fdr': fdr,
        "ge_enrichment_rank_direction": -1,
        "caller_identity": caller_identity,
    }
    
    try:
        # Write the configuration file
        with open(config_path, 'wb') as f:
            tomli_w.dump(config, f)
        
        logger.info(f"Created initial configuration file at {config_path}")
        return config_path
    except Exception as e:
        logger.error(f"Error writing configuration to {config_path}: {e}")
        raise