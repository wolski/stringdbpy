# pyright: reportUnknownMemberType=false
"""FASTA and local-taxonomy species resolvers."""

from __future__ import annotations

import re
import zipfile
from collections import Counter
from importlib import resources
from io import BytesIO
from pathlib import Path
from typing import Protocol

import polars as pl


class SupportedTaxonomy(Protocol):
    """Capability required by FASTA resolvers."""

    def nearest_supported(self, taxon: int) -> int | None:
        """Return the nearest STRING-supported taxon."""
        ...


class LocalTaxonomy:
    """Resolve NCBI taxa to the nearest STRING-supported ancestor."""

    def __init__(self) -> None:
        mapping_root = resources.files("string_gsea.data.mappings")
        with resources.as_file(mapping_root.joinpath("species.v12.0.zip")) as species_path:
            with zipfile.ZipFile(species_path) as zipped:
                species_data = zipped.read("species.v12.0.txt")
            self._string_species = pl.read_csv(BytesIO(species_data), separator="\t").rename(
                {"#taxon_id": "taxon_id"}
            )
        with resources.as_file(mapping_root.joinpath("NCBI_nodes.zip")) as nodes_path:
            with zipfile.ZipFile(nodes_path) as zipped:
                nodes_data = zipped.read("nodes.tsv")
            self._nodes = pl.read_csv(BytesIO(nodes_data), separator="\t")

    def nearest_supported(self, taxon: int) -> int | None:
        supported = set(self._string_species.get_column("taxon_id").to_list())
        current = taxon
        while current != 1:
            if current in supported:
                return current
            parent = self._nodes.filter(pl.col("taxon_id") == current)
            if parent.is_empty():
                return None
            current = int(parent.get_column("parent_taxon_id")[0])
        return None


def _taxon_from_lines(lines: list[str]) -> int | None:
    pattern = re.compile(r"OX=(\d+)")
    taxa = [
        int(match.group(1))
        for line in lines
        if line.startswith(">") and (match := pattern.search(line)) is not None
    ]
    return Counter(taxa).most_common(1)[0][0] if taxa else None


class FastaArchiveResolver:
    """Resolve species from FASTA entries embedded in a ZIP archive."""

    def __init__(self, archive: Path, taxonomy: SupportedTaxonomy) -> None:
        self._archive = archive
        self._taxonomy = taxonomy

    def resolve(self) -> int | None:
        with zipfile.ZipFile(self._archive) as zipped:
            fasta_files = [name for name in zipped.namelist() if name.endswith((".fas", ".fasta"))]
            lines = [
                line for name in fasta_files for line in zipped.read(name).decode().splitlines()
            ]
        taxon = _taxon_from_lines(lines)
        return None if taxon is None else self._taxonomy.nearest_supported(taxon)


class FastaFileResolver:
    """Resolve species from one FASTA file."""

    def __init__(self, path: Path, taxonomy: SupportedTaxonomy) -> None:
        self._path = path
        self._taxonomy = taxonomy

    def resolve(self) -> int | None:
        taxon = _taxon_from_lines(self._path.read_text().splitlines())
        return None if taxon is None else self._taxonomy.nearest_supported(taxon)
