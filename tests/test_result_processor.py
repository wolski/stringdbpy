"""Tests for write_gsea_xlsx() — GSEA result post-processing."""

from pathlib import Path

from string_gsea.gsea_result_processor import write_gsea_xlsx
from string_gsea.models.gsea_models import parse_gsea_tsv_dir

# Path to real TSV fixture data
FIXTURE_GSEA_DIR = Path(__file__).parent / "data" / "outputs" / "mouse_xlsx" / "WU_mouse_fasta_GSEA"


def test_write_gsea_xlsx_creates_xlsx_files(tmp_path):
    """write_gsea_xlsx should create long, pivoted, and merged XLSX files."""
    if not FIXTURE_GSEA_DIR.exists():
        return

    # Find a subdirectory with TSV files
    subdirs = [p for p in FIXTURE_GSEA_DIR.iterdir() if p.is_dir()]
    if not subdirs:
        return

    gsea_result = parse_gsea_tsv_dir(subdirs[0])
    write_gsea_xlsx(gsea_result, "test123", tmp_path)

    xlsx_files = list(tmp_path.glob("*.xlsx"))
    assert len(xlsx_files) >= 3, f"Expected at least 3 XLSX files (long, pivoted, merged), got {len(xlsx_files)}"

    names = [f.name for f in xlsx_files]
    assert any("long" in n for n in names), "Missing long XLSX"
    assert any("pivoted" in n for n in names), "Missing pivoted XLSX"
    assert any("merged" in n for n in names), "Missing merged XLSX"
