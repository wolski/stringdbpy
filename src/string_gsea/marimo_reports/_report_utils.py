"""Shared utilities for Marimo report generation."""

import os
from dataclasses import dataclass
from pathlib import Path

import polars as pl


@dataclass
class ReportConfig:
    """Configuration for report generation."""

    data_file: Path
    links_file: Path
    fdr_threshold: float = 0.05
    genes_mapped_threshold: int = 10
    max_terms: int = 100

    @classmethod
    def from_env(cls) -> "ReportConfig":
        """Load configuration from environment variables."""
        data_file = os.environ.get("GSEA_DATA_FILE", "")
        links_file = os.environ.get("GSEA_LINKS_FILE", "")
        fdr_threshold = float(os.environ.get("GSEA_FDR_THRESHOLD", "0.05"))
        genes_mapped_threshold = int(os.environ.get("GSEA_GENES_THRESHOLD", "10"))
        max_terms = int(os.environ.get("GSEA_MAX_TERMS", "100"))

        if not data_file:
            raise ValueError("GSEA_DATA_FILE environment variable not set")
        if not links_file:
            raise ValueError("GSEA_LINKS_FILE environment variable not set")

        return cls(
            data_file=Path(data_file),
            links_file=Path(links_file),
            fdr_threshold=fdr_threshold,
            genes_mapped_threshold=genes_mapped_threshold,
            max_terms=max_terms,
        )

    @classmethod
    def from_args(cls) -> "ReportConfig":
        """Load configuration from command-line arguments.

        Expected usage: python report.py data_file links_file [--fdr 0.05] [--genes 10]
        """
        import argparse

        parser = argparse.ArgumentParser(description="GSEA Report Configuration")
        parser.add_argument("data_file", type=Path, help="Path to XLSX data file")
        parser.add_argument("links_file", type=Path, help="Path to links.txt file")
        parser.add_argument(
            "--fdr", type=float, default=0.05, help="FDR threshold (default: 0.05)"
        )
        parser.add_argument(
            "--genes", type=int, default=10, help="Genes mapped threshold (default: 10)"
        )
        parser.add_argument(
            "--max-terms", type=int, default=100, help="Max terms per category/contrast (default: 100)"
        )

        args = parser.parse_args()

        return cls(
            data_file=args.data_file,
            links_file=args.links_file,
            fdr_threshold=args.fdr,
            genes_mapped_threshold=args.genes,
            max_terms=args.max_terms,
        )


def load_data(config: ReportConfig) -> pl.DataFrame:
    """Load and return the GSEA results DataFrame."""
    return pl.read_excel(config.data_file)


def load_links(config: ReportConfig) -> dict[str, str]:
    """Load and return the links mapping from links.txt."""
    with open(config.links_file) as f:
        return {
            line.split(": ", 1)[0]: line.split(": ", 1)[1].strip()
            for line in f.read().splitlines()
            if line.strip()
        }


def get_workunit_id(config: ReportConfig) -> str:
    """Extract workunit ID from data file name.

    Examples:
        WU_123_string_gsea_results_long.xlsx → WU_123
        WU2958180_string_gsea_results_long.xlsx → WU2958180
    """
    name = config.data_file.stem
    if "_string_gsea_results" in name:
        return name.split("_string_gsea_results")[0]
    return name
