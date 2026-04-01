import shutil
from pathlib import Path

import requests
from loguru import logger

from string_gsea.gsea_session import GSEASession


class StringGSEAResults:
    """
    Downloads STRING-db GSEA results and holds them in memory.

    Call ``download()`` to fetch TSVs and PNGs from STRING.
    Then use ``write_tsv()``, ``write_graphs()``, ``write_links()``
    to persist to disk, or access ``tsv_content`` / ``graph_content``
    directly for in-memory processing.
    """

    def __init__(self, session: GSEASession):
        self.session = session

        if not self.session.res_data:
            raise RuntimeError("Session contains no result data. Ensure you've polled before building results.")

        self.tsv_content: dict[tuple[str, str], str] = {}
        self.graph_content: dict[tuple[str, str], bytes] = {}

    def download(self) -> "StringGSEAResults":
        """Fetch TSV and PNG results from STRING, storing them in memory."""
        for key, data in self.session.res_data.items():
            if data.get("status") != "success":
                continue
            if "download_url" in data:
                resp = requests.get(data["download_url"])
                resp.raise_for_status()
                self.tsv_content[key] = resp.text
                logger.info(f"Downloaded TSV for {key}")
            if "graph_url" in data:
                resp = requests.get(data["graph_url"])
                resp.raise_for_status()
                self.graph_content[key] = resp.content
                logger.info(f"Downloaded graph for {key}")
        return self

    def get_links(self) -> dict:
        """Extract page URLs for all successful runs."""
        links = {}
        for (outer, inner), data in self.session.res_data.items():
            if data.get("status") == "success" and "page_url" in data:
                links.setdefault(outer, {})[inner] = data["page_url"]
        return links

    def write_links(self, out_dir: Path) -> None:
        """Write links.txt for each outer key."""
        for outer, inners in self.get_links().items():
            sub = out_dir / outer
            sub.mkdir(parents=True, exist_ok=True)
            file_path = sub / "links.txt"
            with open(file_path, "w") as f:
                for inner, url in inners.items():
                    f.write(f"{inner}: {url}\n")
            logger.info(f"Wrote links file: {file_path}")

    def write_tsv(self, out_dir: Path) -> None:
        """Write downloaded TSV content to disk."""
        for (outer, inner), content in self.tsv_content.items():
            sub = out_dir / outer
            sub.mkdir(parents=True, exist_ok=True)
            file_path = sub / f"{inner}_results.tsv"
            file_path.write_text(content)
            logger.info(f"Saved TSV: {file_path}")

    def write_graphs(self, out_dir: Path) -> None:
        """Write downloaded graph PNGs to disk."""
        for (outer, inner), content in self.graph_content.items():
            sub = out_dir / outer
            sub.mkdir(parents=True, exist_ok=True)
            file_path = sub / f"{inner}_results.png"
            file_path.write_bytes(content)
            logger.info(f"Saved graph: {file_path}")

    @staticmethod
    def zip_folder(folder_path: Path) -> Path:
        """Create a .zip archive of the given folder."""
        archive = shutil.make_archive(str(folder_path.parent / folder_path.name), "zip", root_dir=str(folder_path))
        return Path(archive)
