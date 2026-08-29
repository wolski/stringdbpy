"""Client-owned capabilities for species resolution."""

from collections.abc import Sequence
from typing import Protocol


class SpeciesNotResolved(ValueError):
    """Raised after all configured species resolvers return no answer."""


class SpeciesResolver(Protocol):
    """One independent source of a STRING-compatible taxon identifier."""

    def resolve(self) -> int | None:
        """Return a supported taxon, or ``None`` when this source has no answer."""
        ...


class TaxonLookup(Protocol):
    """STRING identifier-to-taxon lookup used by the taxonomy component."""

    def lookup_taxa(self, identifiers: Sequence[str]) -> dict[str, int | None]:
        """Resolve identifiers to their optional NCBI taxon IDs."""
        ...
