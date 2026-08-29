"""Cyclopts boundary and composition root for GSEA runs."""

from dataclasses import replace
from pathlib import Path
from typing import Literal

from cyclopts import App

from string_gsea.configuration import get_configuration
from string_gsea.gsea.application import RunGSEA, RunGSEARequest
from string_gsea.gsea.input import AnalysisName, RankSourceRequest, select_rank_source
from string_gsea.gsea.model.session import SessionSettings
from string_gsea.gsea.output import GSEAOutputWriter
from string_gsea.stringdb_adapters import RequestsStringDB
from string_gsea.taxonomy.application import OrderedSpeciesResolver
from string_gsea.taxonomy.fasta import FastaArchiveResolver, LocalTaxonomy
from string_gsea.taxonomy.stringdb import IdentifierSpeciesResolver

app = App()


@app.default()
def string_gsea_run(
    zip_path: str,
    workunit_id: str,
    out_dir: str = ".",
    which: Literal[
        "pep_1", "pep_1_no_imputed", "pep_2", "pep_2_no_imputed", "none"
    ] = "pep_2_no_imputed",
    fdr: float | None = None,
    create_zip: bool = False,
) -> None:
    """Run STRING-DB GSEA while preserving the established CLI contract."""
    archive = Path(zip_path)
    if not archive.exists():
        raise FileNotFoundError(f"Zip file not found: {archive}")
    output = Path(out_dir)
    output.mkdir(exist_ok=True)

    config = get_configuration()
    if fdr is not None:
        config = replace(config, fdr=fdr)
    settings = SessionSettings(
        api_key=config.api_key,
        fdr=config.fdr,
        ge_enrichment_rank_direction=config.ge_enrichment_rank_direction,
        caller_identity=config.caller_identity,
        creation_date=config.creation_date,
        api_base_url=config.api_base_url,
    )
    gateway = RequestsStringDB(settings)
    analysis = None if which == "none" else AnalysisName(which)
    ranks = select_rank_source(RankSourceRequest(archive=archive, analysis=analysis))
    species = OrderedSpeciesResolver(
        (
            FastaArchiveResolver(archive, LocalTaxonomy()),
            IdentifierSpeciesResolver(ranks.first().sample_identifiers(), gateway),
        )
    ).resolve()
    RunGSEA(gateway, GSEAOutputWriter()).execute(
        RunGSEARequest(
            ranks=ranks,
            workunit_id=workunit_id,
            output_directory=output,
            species=species,
            settings=settings,
            create_zip=create_zip,
        )
    )


if __name__ == "__main__":
    app()
