# pyright: reportUnknownMemberType=false
"""Tests for the injected ORA use case and output compatibility."""

import json
from collections.abc import Sequence
from importlib import resources
from pathlib import Path
from typing import cast

import pytest
from jsonschema import Draft202012Validator
from jsonschema.exceptions import ValidationError

from string_gsea.ora.application import RunORA, RunORARequest
from string_gsea.ora.model import ORARecord
from string_gsea.ora.output import ORAOutputWriter


class FakeORAGateway:
    """Deterministic implementation of the ORA gateway protocol."""

    def map_identifiers(self, identifiers: Sequence[str], species: int) -> dict[str, str]:
        return {
            identifier: f"{species}.{identifier}"
            for identifier in identifiers
            if identifier != "missing"
        }

    def enrich(
        self,
        identifiers: Sequence[str],
        background: Sequence[str],
        species: int,
    ) -> tuple[ORARecord, ...]:
        return (
            {
                "category": "Process",
                "term": "GO:1",
                "number_of_genes": len(identifiers),
                "number_of_genes_in_background": len(background),
                "ncbiTaxonId": species,
                "inputGenes": list(identifiers),
                "preferredNames": ["P1"],
                "p_value": 0.001,
                "fdr": 0.01,
                "description": "Example process",
            },
        )

    def network_link(self, identifiers: Sequence[str], species: int) -> str | None:
        return f"https://string.test/{species}/{len(identifiers)}"


def test_run_ora_writes_compatible_outputs(tmp_path: Path) -> None:
    artifacts = RunORA(FakeORAGateway(), ORAOutputWriter()).execute(
        RunORARequest(
            significant=("P1", "missing"),
            background=("P1", "P2"),
            species=9606,
            output_directory=tmp_path,
            workunit_id="W1",
        )
    )

    assert artifacts.result_directory == tmp_path / "ORA_W1"
    assert set(artifacts.files) == {"json", "tsv", "links", "mapping"}
    mapping = json.loads(artifacts.files["mapping"].read_text())
    assert mapping["unmapped_significant"] == ["missing"]
    assert artifacts.files["links"].read_text() == "STRING-DB Network: https://string.test/9606/1\n"
    assert "P1" in artifacts.files["tsv"].read_text()

    raw_schema: object = json.loads(
        resources.files("string_gsea.ora").joinpath("enrichment_results.schema.json").read_text()
    )
    if not isinstance(raw_schema, dict):
        raise TypeError("ORA result schema must be a JSON object")
    schema = cast(dict[str, object], raw_schema)
    Draft202012Validator.check_schema(schema)
    validator = Draft202012Validator(schema)
    enrichment: object = json.loads(artifacts.files["json"].read_text())
    validator.validate(enrichment)

    invalid = cast(list[dict[str, object]], enrichment)
    invalid[0].pop("term")
    with pytest.raises(ValidationError, match="term"):
        validator.validate(invalid)


@pytest.mark.parametrize(
    "significant, background, message",
    [
        ((), ("P1",), "No identifiers in significant"),
        (("P1",), (), "No identifiers in background"),
        (("missing",), ("P1",), "No significant identifiers could be mapped"),
    ],
)
def test_run_ora_rejects_missing_required_data(
    tmp_path: Path,
    significant: tuple[str, ...],
    background: tuple[str, ...],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        RunORA(FakeORAGateway(), ORAOutputWriter()).execute(
            RunORARequest(significant, background, 9606, tmp_path, "W1")
        )
