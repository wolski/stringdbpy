"""Concrete filesystem output for a completed GSEA application run."""

from __future__ import annotations

import shutil
from pathlib import Path

from loguru import logger

from string_gsea.gsea.model.enrichment import GSEAResult
from string_gsea.gsea.model.ranks import RankListCollection
from string_gsea.gsea.model.session import DownloadedJob, GSEAArtifacts, GSEASession
from string_gsea.gsea.result_processing import write_gsea_xlsx
from string_gsea.gsea.session_yaml import SessionYaml


class GSEAOutputWriter:
    """Write the established result directory and optional archive."""

    def write(
        self,
        *,
        session: GSEASession,
        ranks: RankListCollection,
        downloads: tuple[DownloadedJob, ...],
        result: GSEAResult,
        create_zip: bool,
    ) -> GSEAArtifacts:
        result_directory = session.base_path / f"WU_{session.workunit_id}_GSEA"
        result_directory.mkdir(parents=True, exist_ok=True)
        session_path = SessionYaml.dump(session, result_directory / "gsea_session.yml")
        self._write_ranks(ranks, result_directory)
        result_json_path = result.to_json(
            result_directory / f"WU{session.workunit_id}_gsea_result.json"
        )
        write_gsea_xlsx(result, session.workunit_id, result_directory)
        self._write_downloads(downloads, result_directory)
        archive_path = self._archive(result_directory) if create_zip else None
        return GSEAArtifacts(
            result_directory=result_directory,
            session_path=session_path,
            result_json_path=result_json_path,
            archive_path=archive_path,
        )

    @staticmethod
    def _write_ranks(ranks: RankListCollection, output: Path) -> None:
        directory = output / ranks.analysis
        directory.mkdir(parents=True, exist_ok=True)
        for rank_list in ranks:
            path = directory / f"{rank_list.contrast}.rnk"
            path.write_text(rank_list.to_rnk_string())
            logger.info("Wrote rank file: {}", path)

    @staticmethod
    def _write_downloads(downloads: tuple[DownloadedJob, ...], output: Path) -> None:
        links: dict[str, list[tuple[str, str]]] = {}
        for downloaded in downloads:
            analysis, contrast = downloaded.key
            directory = output / analysis
            directory.mkdir(parents=True, exist_ok=True)
            if downloaded.tsv is not None:
                (directory / f"{contrast}_results.tsv").write_text(downloaded.tsv)
            if downloaded.graph is not None:
                (directory / f"{contrast}_results.png").write_bytes(downloaded.graph)
            if downloaded.page_url is not None:
                links.setdefault(analysis, []).append((contrast, downloaded.page_url))
        for analysis, entries in links.items():
            content = "".join(f"{contrast}: {url}\n" for contrast, url in entries)
            (output / analysis / "links.txt").write_text(content)

    @staticmethod
    def _archive(directory: Path) -> Path:
        archive = shutil.make_archive(
            str(directory.parent / directory.name), "zip", root_dir=str(directory)
        )
        return Path(archive)
