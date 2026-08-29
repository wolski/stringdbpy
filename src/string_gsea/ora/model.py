"""Typed ORA application values."""

from dataclasses import dataclass
from pathlib import Path

type JsonValue = str | int | float | bool | list[JsonValue] | dict[str, JsonValue] | None
type ORARecord = dict[str, JsonValue]


@dataclass(frozen=True, slots=True)
class IdentifierMapping:
    """Mapped and unmapped identifiers for one ORA run."""

    significant: dict[str, str]
    background: dict[str, str]
    unmapped_significant: tuple[str, ...]
    unmapped_background: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ORAResult:
    """Complete runtime result before filesystem serialization."""

    records: tuple[ORARecord, ...]
    mapping: IdentifierMapping
    network_link: str | None


@dataclass(frozen=True, slots=True)
class ORAArtifacts:
    """Paths written for a completed ORA run."""

    result_directory: Path
    files: dict[str, Path]
