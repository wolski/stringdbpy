"""Typed ORA application use case."""

from dataclasses import dataclass
from pathlib import Path

from string_gsea.ora.model import IdentifierMapping, ORAArtifacts, ORAResult
from string_gsea.ora.output import ORAOutputWriter
from string_gsea.ora.ports import ORAGateway


@dataclass(frozen=True, slots=True)
class RunORARequest:
    """Complete input to one ORA execution."""

    significant: tuple[str, ...]
    background: tuple[str, ...]
    species: int
    output_directory: Path
    workunit_id: str


class RunORA:
    """Map identifiers, run enrichment, and persist compatible artifacts."""

    def __init__(self, gateway: ORAGateway, output: ORAOutputWriter) -> None:
        self._gateway = gateway
        self._output = output

    def execute(self, request: RunORARequest) -> ORAArtifacts:
        if not request.significant:
            raise ValueError("No identifiers in significant file")
        if not request.background:
            raise ValueError("No identifiers in background file")

        significant_mapping = self._gateway.map_identifiers(request.significant, request.species)
        background_mapping = self._gateway.map_identifiers(request.background, request.species)
        significant_ids = tuple(significant_mapping.values())
        background_ids = tuple(background_mapping.values())
        if not significant_ids:
            raise ValueError("No significant identifiers could be mapped to STRING IDs")
        if not background_ids:
            raise ValueError("No background identifiers could be mapped to STRING IDs")

        result = ORAResult(
            records=self._gateway.enrich(significant_ids, background_ids, request.species),
            mapping=IdentifierMapping(
                significant=significant_mapping,
                background=background_mapping,
                unmapped_significant=tuple(
                    identifier
                    for identifier in request.significant
                    if identifier not in significant_mapping
                ),
                unmapped_background=tuple(
                    identifier
                    for identifier in request.background
                    if identifier not in background_mapping
                ),
            ),
            network_link=self._gateway.network_link(significant_ids, request.species),
        )
        return self._output.write(result, request.output_directory, request.workunit_id)
