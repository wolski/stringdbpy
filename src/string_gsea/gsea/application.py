"""Typed GSEA application use case."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from string_gsea.gsea.model.enrichment import RunMetadata, parse_gsea_results
from string_gsea.gsea.model.ranks import RankListCollection
from string_gsea.gsea.model.session import GSEAArtifacts, GSEASession, SessionSettings
from string_gsea.gsea.output import GSEAOutputWriter
from string_gsea.gsea.ports import GSEAGateway


@dataclass(frozen=True, slots=True)
class RunGSEARequest:
    """Complete input to one GSEA execution."""

    ranks: RankListCollection
    workunit_id: str
    output_directory: Path
    species: int
    settings: SessionSettings
    create_zip: bool = False


class RunGSEA:
    """Submit ranks, await typed results, and write compatible artifacts."""

    def __init__(
        self,
        gateway: GSEAGateway,
        output: GSEAOutputWriter,
        now: Callable[[], datetime] = datetime.now,
    ) -> None:
        self._gateway = gateway
        self._output = output
        self._now = now

    def execute(self, request: RunGSEARequest) -> GSEAArtifacts:
        session = GSEASession(
            current_date=self._now().strftime("%Y-%m-%d %H:%M:%S"),
            workunit_id=request.workunit_id,
            species=request.species,
            settings=request.settings,
            base_path=request.output_directory,
        )
        analysis = request.ranks.analysis
        for rank_list in request.ranks:
            key = analysis, rank_list.contrast
            session.job_ids[key] = self._gateway.submit_ranks(
                rank_list.to_rnk_string(), request.species
            )
        for key, job_id in session.job_ids.items():
            session.completed_jobs[key] = self._gateway.wait(job_id)

        downloads = tuple(
            self._gateway.download(key, job) for key, job in session.completed_jobs.items()
        )
        tsv_content = {
            downloaded.key: downloaded.tsv for downloaded in downloads if downloaded.tsv is not None
        }
        links = {
            f"{downloaded.key[1]}_results.tsv": downloaded.page_url
            for downloaded in downloads
            if downloaded.page_url is not None
        }
        metadata = RunMetadata(
            workunit_id=request.workunit_id,
            species=request.species,
            fdr=request.settings.fdr,
            ge_enrichment_rank_direction=request.settings.ge_enrichment_rank_direction,
            caller_identity=request.settings.caller_identity,
            api_base_url=request.settings.api_base_url,
        )
        result = parse_gsea_results(
            request.ranks,
            tsv_content,
            metadata=metadata,
            links=links,
        )
        return self._output.write(
            session=session,
            ranks=request.ranks,
            downloads=downloads,
            result=result,
            create_zip=request.create_zip,
        )
