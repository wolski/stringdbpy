from cyclopts import App
from pathlib import Path
from loguru import logger

from string_gsea.config import write_initial_configuration

app = App()


@app.default()
def write_config( caller_identity: str, fdr: float = 0.25,):
    """
    Write an empty configuration file with default values.
    
    Args:
        caller_identity: A string identifying the caller (e.g. "www.example.com")
        fdr: False discovery rate threshold for enrichment analysis (default: 0.25)
    The configuration file is located in:
    - Windows: %APPDATA%/string_gsea/config.toml
    - macOS/Linux: ~/.config/string_gsea/config.toml
    """
    try:
        config_path = write_initial_configuration(caller_identity, fdr)
        #logger.info(f"Configuration file created successfully at: {config_path}")
        
    except Exception as e:
        logger.error(f"Failed to create configuration file: {e}")
        raise


if __name__ == '__main__':
    app() 