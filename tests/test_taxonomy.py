"""Tests for injected and ordered species resolvers."""

import zipfile
from collections.abc import Sequence
from pathlib import Path

import pytest

from string_gsea.taxonomy.application import OrderedSpeciesResolver
from string_gsea.taxonomy.fasta import FastaArchiveResolver, FastaFileResolver, LocalTaxonomy
from string_gsea.taxonomy.ports import SpeciesNotResolved
from string_gsea.taxonomy.stringdb import IdentifierSpeciesResolver


class IdentityTaxonomy:
    def nearest_supported(self, taxon: int) -> int | None:
        return taxon if taxon != 999 else None


class FakeTaxonLookup:
    def lookup_taxa(self, identifiers: Sequence[str]) -> dict[str, int | None]:
        return {identifier: 9606 if identifier != "unknown" else None for identifier in identifiers}


def test_fasta_file_resolver_uses_most_common_ox(tmp_path: Path) -> None:
    fasta = tmp_path / "proteome.fasta"
    fasta.write_text(">a OX=9606\nAA\n>b OX=9606\nAA\n>c OX=10090\nAA\n")

    assert FastaFileResolver(fasta, IdentityTaxonomy()).resolve() == 9606


def test_fasta_archive_resolver_returns_none_without_evidence(tmp_path: Path) -> None:
    archive = tmp_path / "input.zip"
    with zipfile.ZipFile(archive, "w") as zipped:
        zipped.writestr("proteome.fasta", ">a no-taxon\nAA\n")

    assert FastaArchiveResolver(archive, IdentityTaxonomy()).resolve() is None


def test_packaged_taxonomy_resolves_supported_species() -> None:
    assert LocalTaxonomy().nearest_supported(9606) == 9606


def test_identifier_resolver_uses_majority_taxon() -> None:
    resolver = IdentifierSpeciesResolver(("P1", "P2", "unknown"), FakeTaxonLookup())

    assert resolver.resolve() == 9606


class NoSpecies:
    def resolve(self) -> int | None:
        return None


class MouseSpecies:
    def resolve(self) -> int | None:
        return 10090


def test_ordered_resolver_uses_first_answer_and_reports_exhaustion() -> None:
    assert OrderedSpeciesResolver((NoSpecies(), MouseSpecies())).resolve() == 10090
    with pytest.raises(SpeciesNotResolved):
        OrderedSpeciesResolver((NoSpecies(),)).resolve()
