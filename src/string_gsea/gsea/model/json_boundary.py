"""Runtime shape checks for persisted GSEA JSON documents."""

from collections.abc import Mapping
from typing import cast


class InvalidModelDocument(ValueError):
    """Raised when persisted model JSON does not match its required shape."""


def json_object(value: object, field: str) -> dict[str, object]:
    """Validate and return a string-keyed JSON object."""
    if not isinstance(value, Mapping):
        raise InvalidModelDocument(f"{field} must be an object")
    mapping = cast(Mapping[object, object], value)
    if not all(isinstance(key, str) for key in mapping):
        raise InvalidModelDocument(f"{field} must have string keys")
    return {str(key): item for key, item in mapping.items()}


def json_array(value: object, field: str) -> list[object]:
    """Validate and return a JSON array."""
    if not isinstance(value, list):
        raise InvalidModelDocument(f"{field} must be an array")
    return cast(list[object], value)


def required_string(value: object, field: str) -> str:
    """Validate a required string field."""
    if not isinstance(value, str):
        raise InvalidModelDocument(f"{field} must be a string")
    return value


def required_number(value: object, field: str) -> float:
    """Validate a required numeric field and normalize it to float."""
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise InvalidModelDocument(f"{field} must be numeric")
    return float(value)


def required_integer(value: object, field: str) -> int:
    """Validate a required integer field."""
    if isinstance(value, bool) or not isinstance(value, int):
        raise InvalidModelDocument(f"{field} must be an integer")
    return value


def string_map(value: object, field: str) -> dict[str, str]:
    """Validate an object whose values are strings."""
    document = json_object(value, field)
    return {key: required_string(item, f"{field}.{key}") for key, item in document.items()}
