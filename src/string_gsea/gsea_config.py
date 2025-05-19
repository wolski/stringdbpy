import os
import tomli
import tomli_w
import requests
import json
from pathlib import Path
from loguru import logger
from typing import Tuple, Optional
from datetime import datetime
from dataclasses import dataclass, asdict

from string_gsea.config import GSEAConfig


@dataclass
class GSEAConfig:
    api_key: str
    fdr: float
    ge_enrichment_rank_direction: int
    caller_identity: str
    creation_date: Optional[str] = None

    required = ['api_key', 'fdr', 'ge_enrichment_rank_direction', 'caller_identity']

    @classmethod
    def read_toml(cls, path: Path) -> "GSEAConfig":
        """
        Read TOML from `path`, validate required fields, and return a GSEAConfig instance.
        """
        with open(path, 'rb') as f:
            data = tomli.load(f)
        missing = [key for key in cls.required if key not in data]
        if missing:
            raise ValueError(
                f"Configuration file is missing required keys: {', '.join(missing)}"
            )

        return cls(
            api_key=data['api_key'],
            fdr=data['fdr'],
            ge_enrichment_rank_direction=data['ge_enrichment_rank_direction'],
            caller_identity=data['caller_identity'],
            creation_date=data.get('creation_date')
        )

    @classmethod
    def from_dict(cls, data: dict) -> GSEAConfig:
        """
        Initialize GSEAConfig from dict.
        """
        missing = [key for key in cls.required if key not in data]
        if missing:
            raise ValueError(
                f"Configuration file is missing required keys: {', '.join(missing)}"
            )
        return cls(
            api_key=data['api_key'],
            fdr=data['fdr'],
            ge_enrichment_rank_direction=data['ge_enrichment_rank_direction'],
            caller_identity=data['caller_identity'],
            creation_date=data.get('creation_date')
        )

    def write_toml(self, path: Path) -> None:
        """
        Write the current configuration to TOML at `path`.
        """
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, 'wb') as f:
            tomli_w.dump(asdict(self), f)


def _get_config_path() -> Path:
    """
    Determine the platform-specific config.toml path for string_gsea.
    """
    if os.name == 'nt':  # Windows
        config_dir = Path(os.environ.get('APPDATA', '')) / 'string_gsea'
    else:
        config_dir = Path.home() / '.config' / 'string_gsea'
    return config_dir / 'config.toml'


def _fetch_api_key(url: str = "https://version-12-0.string-db.org/api/json/get_api_key") -> Tuple[str, str]:
    """
    Fetch the API key (and optional note) from STRING-DB API.

    Returns:
        api_key: The fetched API key.
        note: Any note associated with the key retrieval.

    Raises:
        requests.RequestException: On network or HTTP errors.
        ValueError: On missing or invalid API key in the response.
    """
    logger.info("Fetching API key from STRING-DB...")
    response = requests.get(url)
    response.raise_for_status()
    data = response.json()

    if not data or not isinstance(data, list) or len(data) == 0:
        raise ValueError("Invalid response format from STRING-DB API")

    api_key = data[0].get("api_key")
    if not api_key:
        raise ValueError("No API key found in the response")

    note = data[0].get("note", "No note provided")
    return api_key, note


def get_configuration() -> GSEAConfig:
    """
    Get the configuration from a TOML file in a platform-independent way.

    The configuration file is located in:
    - Windows: %APPDATA%/string_gsea/config.toml
    - macOS/Linux: ~/.config/string_gsea/config.toml
    """
    config_path = _get_config_path()
    if not config_path.exists():
        raise FileNotFoundError(
            f"Configuration file not found at {config_path}. "
            f"Please create a configuration file using write_initial_configuration() and update it with your settings."
        )

    # Let any tomli or I/O errors propagate
    return GSEAConfig.read_toml(config_path)


def write_initial_configuration(caller_identity: str = "www.fgcz.ch", fdr: float = 0.25) -> Path:
    """
    Write an initial configuration file with the provided parameters and fetch an API key from STRING-DB.

    The configuration file is located in:
    - Windows: %APPDATA%/string_gsea/config.toml
    - macOS/Linux: ~/.config/string_gsea/config.toml
    """
    config_path = _get_config_path()
    config_path.parent.mkdir(parents=True, exist_ok=True)

    if config_path.exists():
        logger.info(f"Configuration file already exists at {config_path}. Overwrite? (y/n)")
        if input().lower() != 'y':
            logger.info("Exiting without overwriting the existing configuration file.")
            return config_path

    # Fetch and propagate any errors
    api_key, note = _fetch_api_key()
    logger.info(f"Successfully obtained API key: {api_key}")
    logger.info(f"Note: {note}")

    config = GSEAConfig(
        api_key=api_key,
        fdr=fdr,
        ge_enrichment_rank_direction=1,
        caller_identity=caller_identity,
        creation_date=datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    )
    # Let write errors propagate
    config.write_toml(config_path)
    logger.info(f"Created initial configuration file at {config_path}")
    return config_path
