"""Client-owned STRING-DB capability for the ORA application."""

from collections.abc import Sequence
from typing import Protocol

from string_gsea.ora.model import ORARecord


class ORAGatewayError(RuntimeError):
    """Raised when the injected ORA gateway cannot satisfy a request."""


class ORAGateway(Protocol):
    """Smallest external capability exercised by the ORA use case."""

    def map_identifiers(self, identifiers: Sequence[str], species: int) -> dict[str, str]:
        """Map submitted identifiers to STRING identifiers."""
        ...

    def enrich(
        self,
        identifiers: Sequence[str],
        background: Sequence[str],
        species: int,
    ) -> tuple[ORARecord, ...]:
        """Return enrichment records for mapped identifiers."""
        ...

    def network_link(self, identifiers: Sequence[str], species: int) -> str | None:
        """Return the optional STRING network page URL."""
        ...
