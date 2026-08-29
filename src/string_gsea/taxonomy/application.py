"""Species-resolution policy composed from injected resolvers."""

from dataclasses import dataclass

from string_gsea.taxonomy.ports import SpeciesNotResolved, SpeciesResolver


@dataclass(frozen=True, slots=True)
class OrderedSpeciesResolver:
    """Return the first answer from an explicitly ordered resolver chain."""

    resolvers: tuple[SpeciesResolver, ...]

    def resolve(self) -> int:
        for resolver in self.resolvers:
            taxon = resolver.resolve()
            if taxon is not None:
                return taxon
        raise SpeciesNotResolved("No configured resolver found a STRING-compatible species")
