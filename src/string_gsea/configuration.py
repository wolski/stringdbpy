"""Configuration file boundary for STRING-GSEA."""

import os
import tomllib
from collections.abc import Mapping
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import ClassVar, Protocol

import tomli_w
from loguru import logger

STRING_API_BASE_DEFAULT = "https://version-12-0.string-db.org/api"


class ApiKeyProvider(Protocol):
    """Volatile external capability required to initialize configuration."""

    def fetch(self) -> tuple[str, str]:
        """Return an API key and its provider note."""
        ...


@dataclass(frozen=True, slots=True)
class GSEAConfig:
    api_key: str
    fdr: float
    ge_enrichment_rank_direction: int
    caller_identity: str
    creation_date: str | None = None
    api_base_url: str = STRING_API_BASE_DEFAULT

    required: ClassVar[tuple[str, ...]] = (
        "api_key",
        "fdr",
        "ge_enrichment_rank_direction",
        "caller_identity",
    )

    @classmethod
    def _validate(cls, data: Mapping[str, object]) -> None:
        """Validate that all required keys are present."""
        missing = [key for key in cls.required if key not in data]
        if missing:
            raise ValueError(f"Configuration is missing required keys: {', '.join(missing)}")

    @classmethod
    def _from_data(cls, data: Mapping[str, object]) -> "GSEAConfig":
        """Construct a GSEAConfig after validating the TOML boundary types."""
        api_key = data["api_key"]
        fdr = data["fdr"]
        rank_direction = data["ge_enrichment_rank_direction"]
        caller_identity = data["caller_identity"]
        creation_date = data.get("creation_date")
        api_base_url = data.get("api_base_url", STRING_API_BASE_DEFAULT)
        if not isinstance(api_key, str) or not api_key:
            raise ValueError("Configuration api_key must be a non-empty string")
        if isinstance(fdr, bool) or not isinstance(fdr, int | float):
            raise ValueError("Configuration fdr must be numeric")
        if isinstance(rank_direction, bool) or not isinstance(rank_direction, int):
            raise ValueError("Configuration ge_enrichment_rank_direction must be an integer")
        if not isinstance(caller_identity, str) or not caller_identity:
            raise ValueError("Configuration caller_identity must be a non-empty string")
        if creation_date is not None and not isinstance(creation_date, str):
            raise ValueError("Configuration creation_date must be a string or null")
        if not isinstance(api_base_url, str) or not api_base_url:
            raise ValueError("Configuration api_base_url must be a non-empty string")
        return cls(
            api_key=api_key,
            fdr=float(fdr),
            ge_enrichment_rank_direction=rank_direction,
            caller_identity=caller_identity,
            creation_date=creation_date,
            api_base_url=api_base_url,
        )

    @classmethod
    def read_toml(cls, path: Path) -> "GSEAConfig":
        """Read TOML from `path`, validate required fields, and return a GSEAConfig instance."""
        with open(path, "rb") as f:
            data = tomllib.load(f)
        cls._validate(data)
        return cls._from_data(data)

    @classmethod
    def from_dict(cls, data: Mapping[str, object]) -> "GSEAConfig":
        """Initialize GSEAConfig from dict."""
        cls._validate(data)
        return cls._from_data(data)

    def write_toml(self, path: Path) -> None:
        """
        Write the current configuration to TOML at `path`.
        """
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "wb") as f:
            tomli_w.dump(asdict(self), f)


def _get_config_path() -> Path:
    """
    Determine the platform-specific config.toml path for string_gsea.
    """
    if os.name == "nt":  # Windows
        config_dir = Path(os.environ.get("APPDATA", "")) / "string_gsea"
    else:
        config_dir = Path.home() / ".config" / "string_gsea"
    return config_dir / "config.toml"


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


def write_initial_configuration(
    api_keys: ApiKeyProvider,
    caller_identity: str = "www.fgcz.ch",
    fdr: float = 0.25,
) -> Path:
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
        if input().lower() != "y":
            logger.info("Exiting without overwriting the existing configuration file.")
            return config_path

    api_key, note = api_keys.fetch()
    logger.info(f"Successfully obtained API key: {api_key}")
    logger.info(f"Note: {note}")

    config = GSEAConfig(
        api_key=api_key,
        fdr=fdr,
        ge_enrichment_rank_direction=1,
        caller_identity=caller_identity,
        creation_date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    )
    # Let write errors propagate
    config.write_toml(config_path)
    logger.info(f"Created initial configuration file at {config_path}")
    return config_path
