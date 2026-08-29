"""Cyclopts boundary for STRING-GSEA configuration creation."""

from cyclopts import App
from loguru import logger

from string_gsea.configuration import write_initial_configuration
from string_gsea.stringdb_adapters import RequestsApiKeyProvider

app = App()


@app.default()
def write_config(caller_identity: str, fdr: float = 0.25) -> None:
    """Create the user configuration file."""
    path = write_initial_configuration(RequestsApiKeyProvider(), caller_identity, fdr)
    logger.info("Configuration file created successfully at: {}", path)


if __name__ == "__main__":
    app()
