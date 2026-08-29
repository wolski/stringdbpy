# pyright: reportUnknownMemberType=false
"""Concrete filesystem serialization for ORA results."""

import json
from pathlib import Path

import polars as pl

from string_gsea.ora.model import ORAArtifacts, ORAResult


class ORAOutputWriter:
    """Write the established ORA JSON, TSV, link, and mapping files."""

    def write(self, result: ORAResult, base_directory: Path, workunit_id: str) -> ORAArtifacts:
        directory = base_directory / f"ORA_{workunit_id}"
        directory.mkdir(parents=True, exist_ok=True)
        files: dict[str, Path] = {}

        json_path = directory / "enrichment_results.json"
        json_path.write_text(json.dumps(result.records, indent=2))
        files["json"] = json_path

        if result.records:
            dataframe = pl.DataFrame(result.records)
            for column in dataframe.columns:
                if dataframe[column].dtype == pl.List:
                    dataframe = dataframe.with_columns(pl.col(column).list.join(",").alias(column))
            tsv_path = directory / "enrichment_results.tsv"
            dataframe.write_csv(tsv_path, separator="\t")
            files["tsv"] = tsv_path

        if result.network_link is not None:
            links_path = directory / "links.txt"
            links_path.write_text(f"STRING-DB Network: {result.network_link}\n")
            files["links"] = links_path

        mapping_path = directory / "id_mapping.json"
        mapping_path.write_text(
            json.dumps(
                {
                    "significant": result.mapping.significant,
                    "background": result.mapping.background,
                    "unmapped_significant": result.mapping.unmapped_significant,
                    "unmapped_background": result.mapping.unmapped_background,
                },
                indent=2,
            )
        )
        files["mapping"] = mapping_path
        return ORAArtifacts(result_directory=directory, files=files)
