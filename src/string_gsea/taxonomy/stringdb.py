"""Identifier-based species resolution through an injected STRING lookup."""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence

from string_gsea.taxonomy.ports import TaxonLookup


class IdentifierSpeciesResolver:
    """Resolve the majority taxon reported for sampled identifiers."""

    def __init__(self, identifiers: Sequence[str], lookup: TaxonLookup) -> None:
        self._identifiers = tuple(identifiers)
        self._lookup = lookup

    def resolve(self) -> int | None:
        taxa = [
            taxon
            for taxon in self._lookup.lookup_taxa(self._identifiers).values()
            if taxon is not None
        ]
        return Counter(taxa).most_common(1)[0][0] if taxa else None
