"""Cyclopts boundary and composition root for ORA runs."""

from pathlib import Path

from cyclopts import App

from string_gsea.configuration import get_configuration
from string_gsea.gsea.model.session import SessionSettings
from string_gsea.ora.application import RunORA, RunORARequest
from string_gsea.ora.output import ORAOutputWriter
from string_gsea.stringdb_adapters import RequestsStringDB
from string_gsea.taxonomy.application import OrderedSpeciesResolver
from string_gsea.taxonomy.fasta import FastaFileResolver, LocalTaxonomy

app = App()


def _read_identifiers(path: Path) -> tuple[str, ...]:
    return tuple(line.strip() for line in path.read_text().splitlines() if line.strip())


@app.default()
def string_ora_run(
    significant: str,
    background: str,
    fasta: str,
    out_dir: str = ".",
    workunit_id: str = "ORA",
) -> dict[str, Path]:
    """Run STRING ORA while preserving the established CLI contract."""
    significant_path = Path(significant)
    background_path = Path(background)
    fasta_path = Path(fasta)
    for path, label in (
        (significant_path, "Significant"),
        (background_path, "Background"),
        (fasta_path, "FASTA"),
    ):
        if not path.exists():
            raise FileNotFoundError(f"{label} file not found: {path}")

    config = get_configuration()
    settings = SessionSettings(
        api_key=config.api_key,
        fdr=config.fdr,
        ge_enrichment_rank_direction=config.ge_enrichment_rank_direction,
        caller_identity=config.caller_identity,
        creation_date=config.creation_date,
        api_base_url=config.api_base_url,
    )
    species = OrderedSpeciesResolver((FastaFileResolver(fasta_path, LocalTaxonomy()),)).resolve()
    artifacts = RunORA(RequestsStringDB(settings), ORAOutputWriter()).execute(
        RunORARequest(
            significant=_read_identifiers(significant_path),
            background=_read_identifiers(background_path),
            species=species,
            output_directory=Path(out_dir),
            workunit_id=workunit_id,
        )
    )
    return artifacts.files


if __name__ == "__main__":
    app()
